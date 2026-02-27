"""
robustness.py
─────────────
Sensitivity and robustness check utilities.

Public functions
────────────────
run_placebo_test       -- Placebo DiD with fake treatment year
run_leave_one_out      -- TWFE re-estimated dropping each state
run_window_sensitivity -- TWFE across different time windows
run_permutation_test   -- Randomization inference (permuted treatment)
"""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from linearmodels.panel import PanelOLS


# ── Placebo test ──────────────────────────────────────────────────────────────

def run_placebo_test(
    panel: pd.DataFrame,
    outcome: str,
    fake_year: int = 2012,
) -> dict:
    """
    Placebo DiD: restrict to pre-2014 data, assign a fake expansion year.

    Expected result under parallel trends: insignificant estimate.
    Uses only 2014 early-adopters vs never-expanded states.
    """
    df = panel[
        (panel['year'] < 2014) &
        (panel['expansion_year'].isin([2014, 0])) &
        panel[outcome].notna()
    ].copy()

    df['fake_post'] = (df['year'] >= fake_year).astype(int)
    df['treat']     = (df['expansion_year'] == 2014).astype(int)
    df['fake_did']  = df['treat'] * df['fake_post']

    model = smf.ols(
        f'{outcome} ~ treat + fake_post + fake_did', data=df
    ).fit(cov_type='cluster', cov_kwds={'groups': df['state_fips']})

    return {
        'fake_year': fake_year,
        'coef':     model.params['fake_did'],
        'se':       model.bse['fake_did'],
        'pvalue':   model.pvalues['fake_did'],
        'ci_lower': model.conf_int().loc['fake_did', 0],
        'ci_upper': model.conf_int().loc['fake_did', 1],
    }


# ── Leave-one-out ─────────────────────────────────────────────────────────────

def run_leave_one_out(panel: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """
    Re-estimate TWFE dropping one state at a time.

    Returns a DataFrame with columns:
        dropped_state, coef, se, pvalue, deviation
    where deviation = |coef - full_sample_coef|.
    """
    # Full-sample estimate for reference
    df_full = panel[panel[outcome].notna()].copy().set_index(['state_fips', 'year'])
    full_model = PanelOLS(
        df_full[outcome], df_full[['post_expansion']],
        entity_effects=True, time_effects=True,
    ).fit(cov_type='clustered', cluster_entity=True)
    full_coef = full_model.params['post_expansion']

    records = []
    for state in panel['state'].unique():
        df = (
            panel[(panel['state'] != state) & panel[outcome].notna()]
            .copy()
            .set_index(['state_fips', 'year'])
        )
        try:
            m = PanelOLS(
                df[outcome], df[['post_expansion']],
                entity_effects=True, time_effects=True,
            ).fit(cov_type='clustered', cluster_entity=True)
            records.append({
                'dropped_state': state,
                'coef':   m.params['post_expansion'],
                'se':     m.std_errors['post_expansion'],
                'pvalue': m.pvalues['post_expansion'],
            })
        except Exception:
            pass

    df_out = pd.DataFrame(records)
    df_out['deviation'] = (df_out['coef'] - full_coef).abs()
    df_out['full_sample_coef'] = full_coef
    return df_out.sort_values('deviation', ascending=False).reset_index(drop=True)


# ── Window sensitivity ────────────────────────────────────────────────────────

def run_window_sensitivity(
    panel: pd.DataFrame,
    outcome: str,
    windows: list = None,
) -> pd.DataFrame:
    """
    Re-estimate TWFE over different time windows.

    Parameters
    ----------
    windows : list of (start_year, end_year) tuples.
              Defaults to six pre-specified windows.
    """
    if windows is None:
        windows = [
            (2010, 2017),   # Short post
            (2010, 2019),   # Medium post
            (2010, 2022),   # Full sample
            (2011, 2022),   # Shorter pre
            (2012, 2022),   # Even shorter pre
            (2010, 2020),   # Exclude COVID
        ]

    records = []
    for start, end in windows:
        df = (
            panel[panel['year'].between(start, end) & panel[outcome].notna()]
            .copy()
            .set_index(['state_fips', 'year'])
        )
        try:
            m = PanelOLS(
                df[outcome], df[['post_expansion']],
                entity_effects=True, time_effects=True,
            ).fit(cov_type='clustered', cluster_entity=True)
            records.append({
                'window':   f'{start}–{end}',
                'coef':     m.params['post_expansion'],
                'se':       m.std_errors['post_expansion'],
                'ci_lower': m.conf_int().loc['post_expansion', 'lower'],
                'ci_upper': m.conf_int().loc['post_expansion', 'upper'],
                'pvalue':   m.pvalues['post_expansion'],
                'n_obs':    m.nobs,
            })
        except Exception as exc:
            print(f"Window {start}–{end} failed: {exc}")

    return pd.DataFrame(records)


# ── Permutation / randomization inference ────────────────────────────────────

def run_permutation_test(
    panel: pd.DataFrame,
    outcome: str,
    n_permutations: int = 500,
    seed: int = 42,
) -> tuple:
    """
    Randomization inference: randomly reassign treatment and re-estimate DiD.

    Returns
    -------
    (actual_estimate, placebo_distribution, p_value)
        placebo_distribution : np.ndarray of length n_permutations
        p_value : proportion of |placebo| >= |actual|
    """
    rng = np.random.default_rng(seed)

    # Actual estimate (simple OLS, not TWFE, for speed)
    df_actual = panel[
        panel['expansion_year'].isin([2014, 0]) &
        panel[outcome].notna()
    ].copy()
    df_actual['post']  = (df_actual['year'] >= 2014).astype(int)
    df_actual['treat'] = (df_actual['expansion_year'] == 2014).astype(int)
    df_actual['did']   = df_actual['treat'] * df_actual['post']
    actual_model   = smf.ols(f'{outcome} ~ treat + post + did', data=df_actual).fit()
    actual_estimate = actual_model.params['did']

    # Null distribution
    state_info = panel.groupby('state_fips')['ever_expanded'].first()
    n_treated  = int(state_info.sum())
    all_states = state_info.index.tolist()

    placebo_dist = []
    for _ in range(n_permutations):
        fake_treated = rng.choice(all_states, size=n_treated, replace=False)
        df = panel[panel[outcome].notna()].copy()
        df['fake_treat'] = df['state_fips'].isin(fake_treated).astype(int)
        df['fake_post']  = (df['year'] >= 2014).astype(int)
        df['fake_did']   = df['fake_treat'] * df['fake_post']
        try:
            m   = smf.ols(f'{outcome} ~ fake_treat + fake_post + fake_did', data=df).fit()
            placebo_dist.append(m.params['fake_did'])
        except Exception:
            pass

    placebo_dist = np.array(placebo_dist)
    p_value = (np.abs(placebo_dist) >= np.abs(actual_estimate)).mean()

    return actual_estimate, placebo_dist, p_value
