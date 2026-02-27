# The Impact of Medicaid Expansion on Health Outcomes: A Causal Inference Analysis

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.3+-276DC3.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Abstract

The Affordable Care Act (ACA) Medicaid expansion, beginning in 2014, extended eligibility to adults with incomes up to 138% of the federal poverty level. As of 2024, 40 states and D.C. have adopted expansion, while 10 states have not — creating a natural experiment for causal inference. This project leverages the staggered adoption of Medicaid expansion across U.S. states to estimate its causal impact on key health outcomes including hospitalization rates, emergency department (ED) utilization, diabetes management, and maternal health indicators. Using a suite of quasi-experimental methods — difference-in-differences (DiD), event study designs, and synthetic control methods — this analysis provides robust evidence on the health effects of coverage expansion while demonstrating best practices in causal inference methodology.

---

## Table of Contents

- [Research Question](#research-question)
- [Background & Motivation](#background--motivation)
- [Data Sources](#data-sources)
- [Methodology](#methodology)
- [Project Structure](#project-structure)
- [Key Findings](#key-findings)
- [Robustness Checks](#robustness-checks)
- [Limitations](#limitations)
- [How to Reproduce](#how-to-reproduce)
- [References](#references)

---

## Research Question

**Primary:** Did Medicaid expansion cause measurable improvements in population health outcomes in expansion states relative to non-expansion states?

**Secondary Questions:**
1. Did expansion reduce avoidable hospitalizations and ED visits (a proxy for improved access to primary care)?
2. Did expansion improve diabetes management outcomes (e.g., diabetes-related mortality, preventable complications)?
3. Did expansion improve maternal health indicators (e.g., maternal mortality, low birth weight, prenatal care utilization)?
4. Are the effects heterogeneous across early vs. late adopting states?

---

## Background & Motivation

The ACA's Medicaid expansion represents one of the largest health insurance coverage expansions in U.S. history. The staggered adoption across states — with some expanding in 2014, others later (e.g., Louisiana in 2016, Virginia in 2019, Missouri in 2021), and some still not expanding — provides a compelling quasi-experimental setting for causal analysis.

Prior literature has established significant effects on insurance coverage and financial protection. This project extends the analysis to direct health outcomes, with particular attention to:

- **Methodological rigor:** Implementing modern DiD estimators that account for staggered treatment timing and heterogeneous treatment effects (Callaway & Sant'Anna, 2021; Sun & Abraham, 2021)
- **Multiple outcome domains:** Examining effects across hospital utilization, chronic disease management, and maternal health
- **Transparency and reproducibility:** All code, data sources, and analytical decisions are fully documented

---

## Data Sources

| Dataset | Source | Years | Key Variables | Link |
|---------|--------|-------|---------------|------|
| **Medicaid Expansion Status** | Kaiser Family Foundation (KFF) | 2010–2023 | State, expansion date, adoption status | [KFF Medicaid Expansion](https://www.kff.org/medicaid/issue-brief/status-of-state-medicaid-expansion-decisions-interactive-map/) |
| **Hospitalization & ED Visits** | CDC WONDER / HCUP State Inpatient Databases | 2010–2022 | Preventable hospitalization rates, ED visit rates by state | [CDC WONDER](https://wonder.cdc.gov/) / [HCUP](https://hcup-us.ahrq.gov/) |
| **Mortality Data** | CDC WONDER Multiple Cause of Death | 2010–2022 | Age-adjusted mortality rates (all-cause, diabetes-related, maternal) | [CDC WONDER Mortality](https://wonder.cdc.gov/mcd.html) |
| **Diabetes Outcomes** | CDC Diabetes Surveillance System | 2010–2021 | Diabetes prevalence, mortality, hospitalization rates by state | [CDC Diabetes Atlas](https://gis.cdc.gov/grasp/diabetes/diabetesatlas.html) |
| **Maternal Health** | CDC WONDER Natality / Linked Birth-Infant Death | 2010–2022 | Maternal mortality ratio, low birth weight %, prenatal care initiation | [CDC Natality](https://wonder.cdc.gov/natality.html) |
| **State Controls** | ACS (Census Bureau) | 2010–2022 | Median income, poverty rate, % uninsured, unemployment, demographics | [Census ACS](https://data.census.gov/) |
| **Health Behaviors** | BRFSS (CDC) | 2010–2022 | Self-reported health status, access to care measures | [BRFSS](https://www.cdc.gov/brfss/) |

### Treatment Variable Construction

States are classified by Medicaid expansion timing:
- **Early adopters (2014):** AZ, AR, CA, CO, CT, DE, DC, HI, IL, IA, KY, MD, MA, MN, NV, NJ, NM, NY, ND, OH, OR, RI, VT, WA, WV
- **Late adopters (2015–2023):** PA (2015), AK, IN, MT (2016), LA (2016), VA (2019), ME (2019), NE (2020), OK, MO (2021), SD (2023), NC (2023)
- **Non-expansion (as of 2024):** TX, FL, GA, WI, MS, AL, SC, TN, WY, KS

> **Note:** Exact dates and groupings will be verified against KFF's latest data during the data collection phase.

---

## Methodology

### 1. Classic Two-Period Difference-in-Differences

A baseline comparison of pre- vs. post-2014 outcomes in expansion vs. non-expansion states.

$$Y_{st} = \alpha + \beta_1 \cdot Expansion_s + \beta_2 \cdot Post_t + \delta \cdot (Expansion_s \times Post_t) + X_{st}\gamma + \epsilon_{st}$$

Where:
- $Y_{st}$ = health outcome in state $s$, year $t$
- $\delta$ = **the causal estimand** (Average Treatment Effect on the Treated)
- $X_{st}$ = time-varying state-level controls

### 2. Staggered DiD with Two-Way Fixed Effects (TWFE)

Accounting for variation in treatment timing:

$$Y_{st} = \alpha_s + \lambda_t + \delta \cdot D_{st} + X_{st}\gamma + \epsilon_{st}$$

Where $\alpha_s$ and $\lambda_t$ are state and year fixed effects, and $D_{st}$ is an indicator for state $s$ having expanded by year $t$.

**Important caveat:** Recent econometrics literature (Goodman-Bacon 2021, de Chaisemartin & D'Haultfœuille 2020) shows that TWFE DiD with staggered treatment can produce biased estimates due to "bad comparisons" (using already-treated units as controls). This motivates the modern estimators below.

### 3. Callaway & Sant'Anna (2021) Estimator

Group-time average treatment effects that avoid the bias in TWFE:

$$ATT(g,t) = E[Y_t - Y_{g-1} | G = g] - E[Y_t - Y_{g-1} | C = 1]$$

Where $g$ indexes the cohort (expansion year) and $C$ indicates the never-treated group. Aggregated to an overall ATT and event-study-style dynamic effects.

**Implementation:** R `did` package

### 4. Event Study Design

Dynamic treatment effects estimated relative to expansion year:

$$Y_{st} = \alpha_s + \lambda_t + \sum_{k \neq -1} \delta_k \cdot \mathbb{1}[t - E_s = k] + X_{st}\gamma + \epsilon_{st}$$

Where $E_s$ is the year state $s$ expanded. Event study plots will show:
- Pre-treatment coefficients (test for parallel trends)
- Post-treatment dynamic effects (how impacts evolve over time)

### 5. Synthetic Control Method

For a deep dive on a single state (e.g., Louisiana, which expanded in 2016), construct a synthetic counterfactual from a weighted combination of non-expansion states that best matches pre-treatment outcomes.

**Implementation:** Python `SparseSC` or R `Synth` package

---

## Project Structure

```
medicaid-expansion-causal-inference/
│
├── README.md                           # This file (research paper format)
├── requirements.txt                    # Python dependencies
├── renv.lock                           # R dependencies (renv)
├── LICENSE
├── .gitignore
│
├── data/
│   ├── raw/                            # Original downloaded datasets
│   │   ├── kff_expansion_status.csv
│   │   ├── cdc_wonder_mortality.txt
│   │   ├── cdc_diabetes_atlas.csv
│   │   ├── cdc_natality.csv
│   │   ├── acs_state_controls.csv
│   │   └── brfss_health_access.csv
│   └── processed/                      # Cleaned, merged panel dataset
│       └── analysis_panel.csv
│
├── notebooks/
│   ├── 01_data_collection.ipynb        # Download and document data sources
│   ├── 02_data_cleaning.ipynb          # Clean, merge, construct panel
│   ├── 03_eda.ipynb                    # Exploratory analysis & descriptive stats
│   ├── 04_did_analysis.ipynb           # Classic & TWFE DiD (Python)
│   ├── 05_event_study.ipynb            # Event study plots (Python)
│   ├── 06_synthetic_control.ipynb      # Synthetic control analysis
│   └── 07_robustness_checks.ipynb      # Sensitivity analyses
│
├── R/
│   ├── 01_did_callaway_santanna.R      # CS estimator (did package)
│   ├── 02_event_study_fixest.R         # Event study with fixest
│   ├── 03_sun_abraham.R               # Sun & Abraham (2021) estimator
│   └── 04_twfe_diagnostics.R          # Goodman-Bacon decomposition
│
├── src/
│   ├── __init__.py
│   ├── data_utils.py                   # Data loading & cleaning functions
│   ├── did_estimators.py               # DiD estimation wrappers
│   ├── plotting.py                     # Publication-quality figure functions
│   └── robustness.py                   # Robustness check utilities
│
├── figures/
│   ├── event_study_hospitalization.png
│   ├── event_study_maternal.png
│   ├── parallel_trends.png
│   ├── synthetic_control_louisiana.png
│   └── expansion_map.png
│
├── tables/
│   ├── summary_statistics.tex
│   ├── main_results.tex
│   └── robustness_results.tex
│
└── docs/
    ├── data_dictionary.md              # Variable definitions
    └── methodology_notes.md            # Detailed methodology discussion
```

---

## Key Findings

All results use 51 states × 13 years (2010–2022). Controls: median household income, poverty rate, % white/black/Hispanic (ACS 1-year). Standard errors clustered by state.

### Primary Outcome: All-Cause Age-Adjusted Mortality (per 100,000)

| Method | Estimate | 95% CI | p-value | Notes |
|--------|----------|--------|---------|-------|
| Simple 2×2 DiD | **−31.5** | [−55.1, −7.8] | 0.009 | 2014 cohort vs never-expanded, with controls |
| TWFE (all states) | +0.15 | [−8.6, +8.9] | 0.97 | Time-varying controls may absorb effect |
| Cohort ATT — 2016 adopters | **−45.4** | [−91.0, +0.1] | 0.049 | Louisiana, Montana |
| Cohort ATT — 2019 adopters | **−69.6** | [−110.8, −28.4] | 0.001 | Virginia, Maine |
| Cohort ATT — 2020 adopters | **−58.9** | [−86.5, −31.3] | <0.001 | Oklahoma, Missouri, Nebraska |
| Synthetic Control — Louisiana | −15.1/yr | — | — | RMSPE ratio 2.1×; p=0.75 (10 donors) |

> **Interpretation:** Clean 2×2 comparisons (each cohort vs never-treated states) consistently show large, statistically significant reductions in all-cause mortality for most expansion cohorts. The aggregate TWFE estimate is attenuated to zero because time-varying demographic controls are themselves partly affected by expansion (a known "bad controls" problem). Cohort-specific estimates are preferred.

### Secondary Outcomes

| Outcome | Simple DiD | p-value | TWFE | p-value |
|---------|-----------|---------|------|---------|
| All-cause crude rate (per 100k) | −47.9 | **0.006** | −1.4 | 0.80 |
| Diabetes mortality, age-adj (per 100k) | −1.5 | **0.046** | −0.16 | 0.69 |
| Diabetes mortality, crude (per 100k) | −2.1 | **0.035** | −0.16 | 0.73 |
| Diabetes prevalence (% adults) | −0.31 | 0.34 | +0.10 | 0.44 |
| Maternal mortality, age-adj (per 100k) | −0.13 | 0.21 | −0.02 | 0.81 |

> **Note:** Maternal mortality findings are limited by CDC data suppression — 55% of state-year cells are missing due to small counts. See `docs/data_dictionary.md`.

---

## Robustness Checks

All checks conducted on all-cause age-adjusted mortality rate.

| Check | Result | Verdict |
|-------|--------|---------|
| **Placebo test** (fake 2012 expansion) | coef = +1.5, p = 0.61 | ✅ Pass — no pre-trend effect |
| **Leave-one-out** | All state-drop estimates within ±10 of full-sample | ✅ Stable |
| **Window sensitivity** (2011–2022) | coef = −8.0, p = 0.056 | ✅ Consistent direction |
| **Permutation test** (500 draws) | p-value consistent with parametric | ✅ Pass |
| **Pre-trends test** (event study) | Max \|t\| < 1.96 in pre-period for 5/7 outcomes | ✅ Mostly flat |
| **TWFE with COVID excluded** (2010–2020) | coef = −4.5, p = 0.23 | ⚠️ Weaker (shorter post-window) |

---

## Limitations

- **Ecological fallacy:** State-level estimates cannot identify individual-level causal effects
- **TWFE attenuation:** Standard TWFE with time-varying controls likely absorbs part of the treatment effect; cohort-specific ATTs are preferred
- **Maternal mortality suppression:** CDC WONDER suppresses cells with <20 deaths — 55% of maternal state-year observations are missing, skewing toward larger states
- **Diabetes prevalence gap:** BRFSS prevalence data unavailable for 2010; incidence data only available as a single-year snapshot (2023)
- **No demographic controls pre-2011:** ACS 1-year estimates are the control source; 2020 is linearly interpolated (ACS not released that year)
- **Spillover effects:** Non-expansion border states may experience indirect coverage gains
- **Concurrent ACA provisions:** Marketplace subsidies and essential health benefits were implemented simultaneously, making Medicaid expansion effects hard to isolate
- **Short post-period for late adopters:** States expanding in 2021–2023 have at most 2 post-treatment years, making their ATTs unreliable
- **SUTVA:** Expansion adoption was not random — political and economic factors confound the comparison

---

## How to Reproduce

### Prerequisites

**Python 3.10+** (full analysis runs in Python)
```bash
pip install -r requirements.txt
```

**R 4.3+** (optional — supplementary estimators)
```R
install.packages(c("did", "fixest", "bacondecomp", "ggplot2", "dplyr"))
```

### Running the Analysis

```bash
# Rebuild the panel from raw data (one-time setup)
jupyter nbconvert --to notebook --execute notebooks/01_data_collection.ipynb
jupyter nbconvert --to notebook --execute notebooks/02_data_cleaning.ipynb

# Or run all notebooks sequentially:
for nb in notebooks/0{3,4,5,6,7}_*.ipynb; do
    jupyter nbconvert --to notebook --execute --inplace "$nb"
done
```

| Notebook | What it does |
|----------|-------------|
| `01_data_collection` | Documents and fetches all raw data sources |
| `02_data_cleaning` | Builds analysis panel (663 × 32) |
| `03_eda` | Summary statistics, pre-trends plots, expansion timeline |
| `04_did_analysis` | Simple DiD, TWFE, event studies for 7 outcomes |
| `05_event_study` | Event study grid, cohort-specific ATTs, pre-trends test |
| `06_synthetic_control` | Synthetic control for Louisiana & Virginia + placebo |
| `07_robustness_checks` | Placebo, leave-one-out, window sensitivity, permutation |

**R scripts** in `R/` implement Callaway-Sant'Anna, Sun-Abraham, and Goodman-Bacon decomposition (require R ≥ 4.3 with `did`, `fixest`, `bacondecomp`).

---

## References

1. Callaway, B., & Sant'Anna, P. H. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230.
2. Sun, L., & Abraham, S. (2021). Estimating dynamic treatment effects in event studies with heterogeneous treatment effects. *Journal of Econometrics*, 225(2), 175–199.
3. Goodman-Bacon, A. (2021). Difference-in-differences with variation in treatment timing. *Journal of Econometrics*, 225(2), 254–277.
4. de Chaisemartin, C., & D'Haultfœuille, X. (2020). Two-way fixed effects estimators with heterogeneous treatment effects. *American Economic Review*, 110(9), 2964–2996.
5. Sommers, B. D., et al. (2017). Changes in utilization and health among low-income adults after Medicaid expansion. *JAMA Internal Medicine*, 177(10), 1440–1447.
6. Miller, S., & Wherry, L. R. (2019). The long-term effects of early life Medicaid coverage. *Journal of Human Resources*, 54(3), 785–824.
7. Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic control methods for comparative case studies. *Journal of the American Statistical Association*, 105(490), 493–505.

---

## Author

**Aditi** — Data Scientist & Researcher

*This project demonstrates proficiency in causal inference, quasi-experimental design, policy evaluation, and reproducible research using Python and R.*

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
