"""
plotting.py
───────────
Publication-quality figure functions for the Medicaid expansion analysis.

Public functions
────────────────
plot_event_study       -- Event study coefficient plot with CI band
plot_pre_trends        -- Outcome trends for expansion vs non-expansion states
plot_expansion_timeline -- Horizontal bar chart of state expansion timing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Shared style ──────────────────────────────────────────────────────────────

COLORS = {
    'expansion':     '#2171b5',
    'no_expansion':  '#cb181d',
    'late':          '#238b45',
    'highlight':     '#e6550d',
    'neutral':       '#969696',
    'ci_fill':       '#9ecae1',
}

_BASE_RCPARAMS = {
    'figure.dpi':       150,
    'font.size':        11,
    'axes.spines.top':  False,
    'axes.spines.right': False,
    'axes.grid':        True,
    'grid.alpha':       0.3,
}

plt.rcParams.update(_BASE_RCPARAMS)


# ── Event study plot ──────────────────────────────────────────────────────────

def plot_event_study(
    es_df: pd.DataFrame,
    outcome: str,
    title: str = None,
    figsize: tuple = (12, 6),
) -> plt.Figure:
    """
    Publication-quality event study plot.

    Parameters
    ----------
    es_df   : DataFrame from run_event_study() with columns
              [event_time, coef, ci_lower, ci_upper]
    outcome : variable name used in the axis label
    title   : optional figure title override
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Confidence band
    ax.fill_between(
        es_df['event_time'], es_df['ci_lower'], es_df['ci_upper'],
        alpha=0.2, color=COLORS['expansion'], label='95% CI',
    )

    # Point estimates + line
    ax.plot(
        es_df['event_time'], es_df['coef'], 'o-',
        color=COLORS['expansion'], linewidth=2, markersize=6,
        label='Coefficient',
    )

    # Reference lines
    ax.axhline(y=0, color='black', linewidth=0.8, zorder=1)
    ax.axvline(x=-0.5, color=COLORS['highlight'], linestyle='--',
               alpha=0.7, linewidth=1.5, label='Expansion')

    # Shade pre-treatment region
    x_min = es_df['event_time'].min()
    ax.axvspan(x_min, -0.5, alpha=0.04, color='gray')

    # Annotations
    yrange = es_df['ci_upper'].max() - es_df['ci_lower'].min()
    y_text = es_df['ci_lower'].min() - 0.06 * yrange
    ax.text(-2.5, y_text, 'Pre-treatment', ha='center',
            fontsize=9, color='gray', style='italic')
    ax.text(3, y_text, 'Post-treatment', ha='center',
            fontsize=9, color='gray', style='italic')

    if title is None:
        title = outcome.replace('_', ' ').title()
    ax.set_title(f'Event Study: {title}', fontsize=14, fontweight='bold')
    ax.set_xlabel('Years Relative to Medicaid Expansion', fontsize=12)
    ax.set_ylabel('Estimated Effect (95% CI)', fontsize=12)
    ax.legend(fontsize=10)

    plt.tight_layout()
    return fig


# ── Pre-trends plot ───────────────────────────────────────────────────────────

def plot_pre_trends(
    panel: pd.DataFrame,
    outcome: str,
    title: str = None,
    figsize: tuple = (12, 6),
) -> plt.Figure | None:
    """
    Plot outcome time series for expansion vs non-expansion states
    with 95% confidence bands.

    Shows the parallel-trends assumption visually.
    """
    if outcome not in panel.columns or panel[outcome].notna().sum() < 20:
        print(f"Not enough data to plot pre-trends for {outcome}")
        return None

    fig, ax = plt.subplots(figsize=figsize)

    for exp_status, color, label in [
        (1, COLORS['expansion'],    'Expansion States'),
        (0, COLORS['no_expansion'], 'Non-Expansion States'),
    ]:
        sub = panel[panel['ever_expanded'] == exp_status]
        g = (
            sub.groupby('year')[outcome]
            .agg(['mean', 'std', 'count'])
            .reset_index()
        )
        g['se']       = g['std'] / np.sqrt(g['count'])
        g['ci_lower'] = g['mean'] - 1.96 * g['se']
        g['ci_upper'] = g['mean'] + 1.96 * g['se']

        ax.plot(g['year'], g['mean'], 'o-', color=color,
                linewidth=2.5, label=label, markersize=5)
        ax.fill_between(g['year'], g['ci_lower'], g['ci_upper'],
                        color=color, alpha=0.12)

    # Treatment reference line
    ax.axvline(x=2014, color='black', linestyle='--', alpha=0.6,
               linewidth=1.5, label='ACA Expansion Begins (2014)')
    ymin, ymax = ax.get_ylim()
    ax.text(2014.15, ymax - (ymax - ymin) * 0.05, 'ACA\nExpansion',
            fontsize=9, va='top', ha='left', style='italic', color='gray')

    if title is None:
        title = outcome.replace('_', ' ').title()
    ax.set_title(f'Pre-Trends Analysis: {title}', fontsize=13, fontweight='bold')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel(title, fontsize=12)
    ax.legend(fontsize=10)

    plt.tight_layout()
    return fig


# ── Expansion timeline plot ───────────────────────────────────────────────────

def plot_expansion_timeline(
    panel: pd.DataFrame,
    figsize: tuple = (14, 16),
) -> plt.Figure:
    """
    Horizontal bar chart showing when each state expanded Medicaid.
    Expanded states are ordered by expansion year; non-expanded at bottom.
    """
    state_info = (
        panel.groupby('state')
        .agg(expansion_year=('expansion_year', 'first'),
             ever_expanded=('ever_expanded', 'first'))
        .reset_index()
    )

    expanded     = state_info[state_info['ever_expanded'] == 1].sort_values('expansion_year')
    not_expanded = state_info[state_info['ever_expanded'] == 0].sort_values('state')
    ordered      = pd.concat([expanded, not_expanded]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=figsize)

    for i, (_, row) in enumerate(ordered.iterrows()):
        if row['ever_expanded'] == 1:
            color = (COLORS['expansion'] if row['expansion_year'] == 2014
                     else COLORS['late'])
            width = 2023 - row['expansion_year']
            ax.barh(i, width, left=row['expansion_year'],
                    color=color, alpha=0.75, height=0.7)
            ax.plot(row['expansion_year'], i, 'o',
                    color=color, markersize=5)
        else:
            ax.barh(i, 0.4, left=2023, color=COLORS['no_expansion'],
                    alpha=0.75, height=0.7)

    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels(ordered['state'], fontsize=8)
    ax.set_xlim(2013, 2024)
    ax.set_xlabel('Year', fontsize=11)
    ax.set_title('Medicaid Expansion Timeline by State',
                 fontsize=14, fontweight='bold')
    ax.axvline(x=2014, color='black', linestyle='--',
               alpha=0.5, linewidth=1.2)

    legend_elements = [
        mpatches.Patch(color=COLORS['expansion'],    alpha=0.75, label='Expanded 2014'),
        mpatches.Patch(color=COLORS['late'],         alpha=0.75, label='Expanded after 2014'),
        mpatches.Patch(color=COLORS['no_expansion'], alpha=0.75, label='Never expanded'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

    plt.tight_layout()
    return fig
