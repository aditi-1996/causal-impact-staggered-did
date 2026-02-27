"""
Medicaid Expansion Causal Inference — Analysis Module
======================================================

Utility functions for data loading, DiD estimation, visualization,
and robustness checks.
"""

from .data_utils import load_panel, parse_cdc_wonder_csv, parse_natality_csv, fetch_acs_data
from .did_estimators import run_simple_did, run_twfe, run_event_study
from .plotting import plot_event_study, plot_pre_trends, plot_expansion_timeline
from .robustness import run_placebo_test, run_leave_one_out, run_window_sensitivity, run_permutation_test

__all__ = [
    "load_panel",
    "parse_cdc_wonder_csv",
    "parse_natality_csv",
    "fetch_acs_data",
    "run_simple_did",
    "run_twfe",
    "run_event_study",
    "plot_event_study",
    "plot_pre_trends",
    "plot_expansion_timeline",
    "run_placebo_test",
    "run_leave_one_out",
    "run_window_sensitivity",
    "run_permutation_test",
]
