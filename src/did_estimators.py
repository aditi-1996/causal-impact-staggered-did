"""
did_estimators.py
─────────────────
Difference-in-differences estimation wrappers.

Public functions
────────────────
run_simple_did   -- Classic 2×2 DiD (early adopters vs never-expanded)
run_twfe         -- Two-way fixed effects DiD (all states, staggered)
run_event_study  -- TWFE event study with leads/lags
"""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from linearmodels.panel import PanelOLS


# ── Simple 2×2 DiD ────────────────────────────────────────────────────────────

def run_simple_did(panel: pd.DataFrame, outcome: str, controls: list = None) -> object:
    """
    Classic 2×2 DiD using early adopters (2014) vs never-expanded states.

    Specification:
        Y = α + β1·Treat + β2·Post + δ·(Treat×Post) + Xγ + ε
        δ is the ATT estimate.

    Standard errors clustered by state.

    Parameters
    ----------
    panel    : analysis panel DataFrame
    outcome  : name of the outcome column
    controls : list of time-varying covariate column names (optional)

    Returns
    -------
    statsmodels RegressionResultsWrapper
    """
    df = panel[
        panel['expansion_year'].isin([2014, 0]) &
        panel[outcome].notna()
    ].copy()

    df['post']        = (df['year'] >= 2014).astype(int)
    df['treat']       = (df['expansion_year'] == 2014).astype(int)
    df['treat_x_post'] = df['treat'] * df['post']

    formula = f'{outcome} ~ treat + post + treat_x_post'
    if controls:
        available = [c for c in controls if c in df.columns and df[c].notna().sum() > 50]
        if available:
            formula += ' + ' + ' + '.join(available)

    model = smf.ols(formula, data=df).fit(
        cov_type='cluster',
        cov_kwds={'groups': df['state_fips']},
    )

    print(f"\n{'='*60}")
    print(f"Simple 2×2 DiD: {outcome}")
    print(f"{'='*60}")
    print(f"Treatment: 2014 expansion | Control: never-expanded")
    print(f"N = {model.nobs:.0f}   R² = {model.rsquared:.4f}")
    coef = model.params['treat_x_post']
    se   = model.bse['treat_x_post']
    ci   = model.conf_int().loc['treat_x_post']
    pv   = model.pvalues['treat_x_post']
    print(f"\n*** DiD Estimate (treat×post): {coef:.4f}  SE: {se:.4f}")
    print(f"    95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]   p = {pv:.4f}")
    return model


# ── TWFE DiD ─────────────────────────────────────────────────────────────────

def run_twfe(panel: pd.DataFrame, outcome: str, controls: list = None) -> object:
    """
    Two-way fixed effects DiD with state and year fixed effects.

    Specification:
        Y_{st} = α_s + λ_t + δ·D_{st} + Xγ + ε_{st}
        D_{st} = post_expansion indicator

    Standard errors clustered by state.

    Note: TWFE can be biased under staggered adoption with heterogeneous
    treatment effects — see R/01_did_callaway_santanna.R for a robust
    alternative.
    """
    df = panel[panel[outcome].notna()].copy()
    df = df.set_index(['state_fips', 'year'])

    y = df[outcome]
    exog_cols = ['post_expansion']
    if controls:
        avail = [c for c in controls if c in df.columns and df[c].notna().sum() > 50]
        exog_cols += avail

    X = df[exog_cols].copy()
    mask = X.notna().all(axis=1) & y.notna()
    y, X = y[mask], X[mask]

    model = PanelOLS(y, X, entity_effects=True, time_effects=True).fit(
        cov_type='clustered', cluster_entity=True
    )

    print(f"\n{'='*60}")
    print(f"TWFE DiD: {outcome}")
    print(f"{'='*60}")
    print(f"All states | Entity + Time FE | N = {model.nobs}")
    print(f"R² (within) = {model.rsquared_within:.4f}")
    coef = model.params['post_expansion']
    se   = model.std_errors['post_expansion']
    ci   = model.conf_int().loc['post_expansion']
    pv   = model.pvalues['post_expansion']
    print(f"\n*** TWFE Estimate (post_expansion): {coef:.4f}  SE: {se:.4f}")
    print(f"    95% CI: [{ci['lower']:.4f}, {ci['upper']:.4f}]   p = {pv:.4f}")
    return model


# ── Event Study ───────────────────────────────────────────────────────────────

def run_event_study(
    panel: pd.DataFrame,
    outcome: str,
    controls: list = None,
    leads: int = 4,
    lags: int = 8,
) -> tuple:
    """
    TWFE event study with leads/lags relative to treatment year.

    Omits t = -1 as the reference period. Endpoints are binned at
    (-leads, +lags) to absorb boundary variation.

    Returns
    -------
    (model, es_df)
        model  : linearmodels PanelOLS result
        es_df  : DataFrame with columns [event_time, coef, se, ci_lower,
                 ci_upper, pvalue] sorted by event_time
    """
    df = panel[
        panel[outcome].notna() &
        (panel['event_time'].notna() | (panel['ever_expanded'] == 0))
    ].copy()

    # Bin endpoints
    df['et_binned'] = df['event_time'].copy()
    df.loc[df['et_binned'] < -leads, 'et_binned'] = -leads
    df.loc[df['et_binned'] >  lags,  'et_binned'] =  lags

    event_times = sorted(
        [t for t in df['et_binned'].dropna().unique() if t != -1]
    )

    for t in event_times:
        col = f'et_{int(abs(t))}m' if t < 0 else f'et_{int(t)}p'
        df[col] = (df['et_binned'] == t).astype(int)

    df_panel = df.set_index(['state_fips', 'year'])
    y = df_panel[outcome]
    et_cols = sorted([c for c in df_panel.columns if c.startswith('et_')])

    exog_cols = et_cols.copy()
    if controls:
        avail = [c for c in controls
                 if c in df_panel.columns and df_panel[c].notna().sum() > 50]
        exog_cols += avail

    X = df_panel[exog_cols].copy()
    mask = X.notna().all(axis=1) & y.notna()
    y, X = y[mask], X[mask]

    model = PanelOLS(y, X, entity_effects=True, time_effects=True).fit(
        cov_type='clustered', cluster_entity=True
    )

    # Extract coefficients
    records = []
    for t in event_times:
        col = f'et_{int(abs(t))}m' if t < 0 else f'et_{int(t)}p'
        if col in model.params.index:
            records.append({
                'event_time': t,
                'coef':       model.params[col],
                'se':         model.std_errors[col],
                'ci_lower':   model.conf_int().loc[col, 'lower'],
                'ci_upper':   model.conf_int().loc[col, 'upper'],
                'pvalue':     model.pvalues[col],
            })

    # Add reference period t = -1 (normalized to 0)
    records.append({'event_time': -1, 'coef': 0, 'se': 0,
                    'ci_lower': 0, 'ci_upper': 0, 'pvalue': np.nan})

    es_df = pd.DataFrame(records).sort_values('event_time').reset_index(drop=True)
    return model, es_df
