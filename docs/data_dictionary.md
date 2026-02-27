# Data Dictionary & Sources Guide

## Overview

This document provides detailed information on each data source used in the analysis, including download instructions, variable definitions, and notes on data quality.

---

## 1. Medicaid Expansion Status (Treatment Variable)

**Source:** Kaiser Family Foundation (KFF)
**URL:** https://www.kff.org/medicaid/issue-brief/status-of-state-medicaid-expansion-decisions-interactive-map/
**Format:** Manual compilation into CSV

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `state` | str | State name |
| `state_fips` | str | 2-digit FIPS code |
| `expansion_status` | str | "Expanded" / "Not Expanded" |
| `expansion_date` | date | Effective date of expansion (YYYY-MM-DD) |
| `expansion_year` | int | Year of expansion (for DiD grouping) |
| `cohort` | str | "Early (2014)" / "Late (2015-2023)" / "Never" |

### Notes
- Some states had partial/early expansions before 2014 (e.g., CT, MN, DC). Decision needed on how to classify these.
- Wisconsin did not formally expand but covers adults up to 100% FPL through a waiver — typically classified as non-expansion.

---

## 2. CDC WONDER — Multiple Cause of Death (Mortality Data)

**Source:** CDC WONDER
**URL:** https://wonder.cdc.gov/mcd.html
**Access:** Web query interface (no API); results exported as tab-delimited .txt

### Query Parameters
- **Group By:** State, Year
- **Years:** 2010–2022
- **Outcomes of interest:**
  - All-cause mortality (age-adjusted rate per 100,000)
  - Diabetes-related mortality (ICD-10: E10-E14)
  - Maternal mortality (ICD-10: O00-O99, A34)

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `state` | str | State of residence |
| `year` | int | Year of death |
| `deaths` | int | Number of deaths |
| `population` | int | State population (denominator) |
| `crude_rate` | float | Deaths per 100,000 |
| `age_adjusted_rate` | float | Age-adjusted rate per 100,000 (using 2000 standard) |

### Notes
- CDC suppresses counts < 10 for privacy. States with small populations may have missing values for specific causes.
- Age-adjusted rates use the 2000 U.S. standard population.

### Known Data Limitation — Maternal Mortality Suppression

**364 of 663 state-year observations (55%) are missing `maternal_age_adj_rate` and `maternal_crude_rate`.** This is an inherent CDC WONDER data quality issue, not a code bug:

- CDC WONDER marks any cell with **fewer than 20 maternal deaths** as `"Unreliable"` and excludes it from downloads, per federal confidentiality standards (45 CFR Part 164).
- The classic WONDER database (2010–2020) covers only **38 of 51 states** for maternal mortality; 13 states have too few deaths in early years to produce reliable rates.
- Most missingness is concentrated in **small-population states** (Alaska, Delaware, Hawaii, Montana, North Dakota, Vermont, Wyoming, etc.) and in **earlier years** (2010–2015) before maternal mortality increased nationally.

**Implication for analysis:** Maternal mortality regressions are restricted to the subset of state-years with non-missing data. This introduces a potential **selection bias** — remaining observations skew toward larger, more populous states. Maternal results should be interpreted with caution and treated as exploratory rather than definitive. As a robustness check, aggregate (pooled) national-level trends can be used to verify directional consistency.

---

## 3. CDC Diabetes Surveillance System

**Source:** CDC Diabetes Atlas
**URL:** https://gis.cdc.gov/grasp/diabetes/diabetesatlas.html
**Format:** Downloadable CSV

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `state` | str | State |
| `year` | int | Year |
| `diabetes_prevalence` | float | % of adults diagnosed with diabetes |
| `diabetes_incidence` | float | New cases per 1,000 adults |
| `diabetes_mortality_rate` | float | Diabetes-related death rate per 100,000 |
| `diabetes_hosp_rate` | float | Diabetes-related hospitalization rate |

### Notes
- Prevalence is based on BRFSS self-report ("Have you ever been told by a doctor that you have diabetes?")
- Hospitalization data may not be available for all states/years

---

## 4. CDC WONDER — Natality (Maternal & Infant Health)

**Source:** CDC WONDER Natality Files
**URL:** https://wonder.cdc.gov/natality.html
**Access:** Web query interface

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `state` | str | State of residence |
| `year` | int | Year of birth |
| `total_births` | int | Total live births |
| `low_birth_weight_pct` | float | % births < 2,500g |
| `preterm_birth_pct` | float | % births < 37 weeks |
| `prenatal_care_first_trimester_pct` | float | % with prenatal care in 1st trimester |
| `maternal_age_mean` | float | Mean maternal age |
| `cesarean_rate` | float | % cesarean deliveries |

### Notes
- Maternal mortality ratios are better sourced from the Maternal Mortality Review data or CDC WONDER mortality files (ICD-10 O-codes) since natality files track births, not deaths.
- The 2003 revised birth certificate was adopted by all states by 2016, creating potential measurement changes in some variables.

---

## 5. American Community Survey (State Controls)

**Source:** U.S. Census Bureau, ACS 1-Year Estimates
**URL:** https://data.census.gov/
**API:** Available via Census API (recommended for reproducibility)

### Variables

| Variable | Type | Description | ACS Table |
|----------|------|-------------|-----------|
| `median_household_income` | float | Median household income ($) | B19013 |
| `poverty_rate` | float | % below federal poverty level | S1701 |
| `uninsured_rate` | float | % without health insurance | S2701 |
| `unemployment_rate` | float | % unemployed (civilian labor force) | S2301 |
| `pct_white` | float | % White non-Hispanic | B03002 |
| `pct_black` | float | % Black non-Hispanic | B03002 |
| `pct_hispanic` | float | % Hispanic/Latino | B03002 |
| `pct_65_plus` | float | % age 65+ | S0101 |
| `pct_female` | float | % female | S0101 |
| `total_population` | int | Total state population | B01001 |

### Notes
- ACS 1-year estimates are available for areas with population ≥ 65,000 (all states qualify).
- Use `census` Python package or direct API calls for reproducible data collection.

---

## 6. BRFSS (Health Behaviors & Access)

**Source:** CDC Behavioral Risk Factor Surveillance System
**URL:** https://www.cdc.gov/brfss/
**Format:** SAS transport files or pre-computed state-level prevalence tables

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `state` | str | State |
| `year` | int | Survey year |
| `has_personal_doctor_pct` | float | % with personal doctor/healthcare provider |
| `could_not_afford_doctor_pct` | float | % who couldn't see doctor due to cost (past 12 months) |
| `routine_checkup_pct` | float | % with routine checkup in past year |
| `self_rated_health_fair_poor_pct` | float | % reporting fair or poor health |

### Notes
- BRFSS is survey-based; estimates have associated confidence intervals.
- State-level prevalence data is available through CDC's interactive data portal.
- Consider using BRFSS SMART data for metro-level analysis if needed.

---

## Panel Dataset Construction

The final analysis panel (`data/processed/analysis_panel.csv`) merges all sources by **state × year** (2010–2022):

```
analysis_panel:
├── Identifiers: state, state_fips, year
├── Treatment: expansion_status, expansion_year, post_expansion, years_since_expansion
├── Outcomes: mortality rates, diabetes indicators, maternal health, ED/hosp rates
├── Controls: income, poverty, uninsured, unemployment, demographics
└── Access: BRFSS health access measures
```

### Key Constructed Variables

| Variable | Definition |
|----------|-----------|
| `post_expansion` | 1 if year ≥ expansion_year, 0 otherwise (0 for never-treated) |
| `years_since_expansion` | year − expansion_year (negative = pre-period; NA for never-treated) |
| `treat` | 1 if state ever expanded by end of study period, 0 otherwise |
| `cohort_group` | Expansion year (2014, 2015, ..., 2023) or 0 for never-treated |
