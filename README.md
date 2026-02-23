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

> *This section will be populated as the analysis progresses.*

| Outcome | Method | Estimate | 95% CI | Significance |
|---------|--------|----------|--------|-------------|
| Preventable Hospitalizations | TWFE DiD | — | — | — |
| ED Visit Rate | Callaway-Sant'Anna | — | — | — |
| Diabetes Mortality | Event Study | — | — | — |
| Maternal Mortality Ratio | Synthetic Control | — | — | — |

---

## Robustness Checks

The following sensitivity analyses will be conducted to validate the main findings:

1. **Parallel Trends Testing** — Formal pre-trend tests and visual inspection of event study pre-period coefficients
2. **Placebo Tests** — Assign fake treatment dates (e.g., 2011) and test for spurious effects
3. **Alternative Control Groups** — Restrict to border-state pairs or demographically similar states
4. **Goodman-Bacon Decomposition** — Decompose TWFE estimate into timing-group sub-estimates to identify potential bias
5. **Leave-One-Out Analysis** — Sequentially drop each state to check for influential observations
6. **Varying Pre/Post Windows** — Test sensitivity to the pre-treatment and post-treatment period definitions
7. **Dose-Response** — Test whether states with larger coverage gains show larger health effects
8. **Covariate Balancing** — Inverse probability weighting to address pre-treatment covariate differences

---

## Limitations

- **Ecological fallacy:** State-level analysis cannot identify individual-level causal effects
- **Data granularity:** Some outcomes are only available annually, limiting precision
- **Spillover effects:** Non-expansion states bordering expansion states may experience indirect effects
- **Concurrent policies:** Other ACA provisions (marketplace subsidies, essential health benefits) were implemented simultaneously
- **Outcome availability:** Some health outcomes have reporting lags or inconsistent state-level coverage
- **SUTVA concerns:** States' expansion decisions are not random — political and economic factors drive adoption

---

## How to Reproduce

### Prerequisites

**Python 3.10+**
```bash
pip install -r requirements.txt
```

**R 4.3+**
```R
# Install required packages
install.packages(c("did", "fixest", "Synth", "bacondecomp", "ggplot2", "dplyr", "modelsummary"))
```

### Running the Analysis

```bash
# Step 1: Data collection & cleaning
jupyter notebook notebooks/01_data_collection.ipynb
jupyter notebook notebooks/02_data_cleaning.ipynb

# Step 2: Exploratory analysis
jupyter notebook notebooks/03_eda.ipynb

# Step 3: Main analysis (Python)
jupyter notebook notebooks/04_did_analysis.ipynb
jupyter notebook notebooks/05_event_study.ipynb
jupyter notebook notebooks/06_synthetic_control.ipynb

# Step 4: Modern DiD estimators (R)
Rscript R/01_did_callaway_santanna.R
Rscript R/02_event_study_fixest.R

# Step 5: Robustness
jupyter notebook notebooks/07_robustness_checks.ipynb
```

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
