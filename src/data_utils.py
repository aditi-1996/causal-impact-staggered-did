"""
data_utils.py
─────────────
Data loading and cleaning utilities for the Medicaid expansion analysis.

Public functions
────────────────
parse_cdc_wonder_csv   -- Combine two-part CDC WONDER CSV exports
parse_natality_csv     -- Parse CDC WONDER Natality CSV
fetch_acs_data         -- Fetch ACS state controls via Census API
load_panel             -- Load the final analysis panel
"""

import os
import requests
import warnings
import numpy as np
import pandas as pd

# ── Column rename maps ────────────────────────────────────────────────────────

_WONDER_RENAME = {
    "State": "state",
    "State Code": "state_fips",
    "Year": "year",
    "Year Code": "year_code",
    "Deaths": "deaths",
    "Population": "population",
    "Crude Rate": "crude_rate",
    "Age Adjusted Rate": "age_adj_rate",
}

_NATALITY_RENAME = {
    "State": "state",
    "State Code": "state_fips",
    "Year": "year",
    "Year Code": "year_code",
    "Births": "total_births",
    "Birth Rate": "birth_rate",
    "Average Birth Weight": "avg_birth_weight_g",
    "% Low Birthweight": "low_birth_weight_pct",
    "% Preterm": "preterm_pct",
    "Total Population": "total_population",
}

_SUPPRESSED = {"Suppressed", "Unreliable", "Not Available", "Missing", "Not Applicable", ""}


# ── CDC WONDER CSV parser ─────────────────────────────────────────────────────

def _read_wonder_csv(filepath: str) -> pd.DataFrame:
    """
    Read a single CDC WONDER CSV export and return clean data rows.

    CDC WONDER CSVs have three types of rows:
      1. Data rows     — Notes column is empty  (the ones we want)
      2. Subtotal rows — Notes column = "Total"
      3. Footer rows   — Notes column contains metadata text

    We keep only rows where Notes is NaN.
    The 2018-2022 expanded database adds confidence-interval columns;
    we drop those to keep a consistent schema.
    """
    df = pd.read_csv(filepath, dtype=str, low_memory=False)

    notes_col = df.columns[0]           # always the first column
    df = df[df[notes_col].isna()].copy()
    df = df.drop(columns=[notes_col])

    # Standardize column names
    df = df.rename(columns={k: v for k, v in _WONDER_RENAME.items() if k in df.columns})

    # Keep only the core columns (drop CI columns present in 2018-2022 files)
    core = ["state", "state_fips", "year", "deaths", "population",
            "crude_rate", "age_adj_rate"]
    df = df[[c for c in core if c in df.columns]].copy()

    # Pad FIPS codes
    df["state_fips"] = df["state_fips"].str.strip().str.zfill(2)

    # Year to int
    df["year"] = pd.to_numeric(df["year"], errors="coerce")

    # Numeric conversions — treat all suppressed/unreliable as NaN
    for col in ["deaths", "population", "crude_rate", "age_adj_rate"]:
        if col in df.columns:
            df[col] = (
                df[col]
                .str.strip()
                .replace(_SUPPRESSED, np.nan)
                .pipe(pd.to_numeric, errors="coerce")
            )

    df = df.dropna(subset=["state_fips", "year"])
    df["year"] = df["year"].astype(int)
    return df


def parse_cdc_wonder_csv(
    file_old: str,
    file_new: str,
    prefix: str,
    old_years: tuple = (2010, 2017),
    new_years: tuple = (2018, 2022),
) -> pd.DataFrame:
    """
    Combine two CDC WONDER CSV exports that split at ~2018.

    CDC WONDER split its databases around 2018: the classic Underlying Cause
    of Death database covers 1999–2020, while the expanded database covers
    2018–present with slightly different population denominators.

    Strategy: use file_old for old_years, file_new for new_years.
    This avoids duplicate rows in the 2018–2020 overlap and keeps each
    year paired with its best population-estimate vintage.

    Parameters
    ----------
    file_old : path to the earlier download (e.g. 2010–2020 file)
    file_new : path to the later download  (e.g. 2018–2022 file)
    prefix   : column prefix, e.g. "allcause", "diabetes", "maternal"
    old_years, new_years : (min, max) year to keep from each file

    Returns
    -------
    DataFrame with columns:
        state_fips, year,
        {prefix}_deaths, {prefix}_population,
        {prefix}_crude_rate, {prefix}_age_adj_rate
    """
    df_old = _read_wonder_csv(file_old)
    df_old = df_old[df_old["year"].between(*old_years)].copy()

    df_new = _read_wonder_csv(file_new)
    df_new = df_new[df_new["year"].between(*new_years)].copy()

    df = pd.concat([df_old, df_new], ignore_index=True)
    df = df.sort_values(["state_fips", "year"]).reset_index(drop=True)

    # Rename to prefixed outcome columns
    rename = {
        "deaths":    f"{prefix}_deaths",
        "population": f"{prefix}_population",
        "crude_rate": f"{prefix}_crude_rate",
        "age_adj_rate": f"{prefix}_age_adj_rate",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df = df.drop(columns=["state"], errors="ignore")

    n_states = df["state_fips"].nunique()
    n_rows   = len(df)
    ymin, ymax = df["year"].min(), df["year"].max()
    pct_missing = df[f"{prefix}_age_adj_rate"].isna().mean() * 100
    print(
        f"  [OK] {prefix:10s}  {n_rows:>4d} rows | {n_states} states "
        f"| years {ymin}-{ymax} | {pct_missing:.0f}% missing age-adj rate"
    )
    return df


# ── Natality parser ───────────────────────────────────────────────────────────

def parse_natality_csv(filepath: str) -> pd.DataFrame:
    """
    Parse CDC WONDER Natality CSV export.

    Keeps data rows only (Notes is NaN). Drops state-level total rows
    (Notes == "Total"). Handles "Not Available" and "Suppressed" values.

    Returns columns: state_fips, year, total_births, birth_rate,
                     avg_birth_weight_g  (and others if present in file)
    """
    df = pd.read_csv(filepath, dtype=str, low_memory=False)

    notes_col = df.columns[0]
    df = df[df[notes_col].isna()].copy()
    df = df.drop(columns=[notes_col])

    df = df.rename(columns={k: v for k, v in _NATALITY_RENAME.items() if k in df.columns})

    df["state_fips"] = df["state_fips"].str.strip().str.zfill(2)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")

    numeric_cols = ["total_births", "birth_rate", "avg_birth_weight_g",
                    "low_birth_weight_pct", "preterm_pct", "total_population"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = (
                df[col]
                .str.strip()
                .replace(_SUPPRESSED, np.nan)
                .pipe(pd.to_numeric, errors="coerce")
            )

    df = df.dropna(subset=["state_fips", "year"])
    df["year"] = df["year"].astype(int)

    # Exclude total_population — it duplicates the ACS control already in the panel
    keep = ["state_fips", "year"] + [
        c for c in numeric_cols if c in df.columns and c != "total_population"
    ]
    df = df[keep].copy()

    n_states = df["state_fips"].nunique()
    print(
        f"  [OK] natality   {len(df):>4d} rows | {n_states} states "
        f"| years {df['year'].min()}-{df['year'].max()}"
    )
    return df


# ── ACS state controls ────────────────────────────────────────────────────────

_ACS_VARS = {
    "B19013_001E": "median_household_income",
    "B17001_001E": "poverty_universe",
    "B17001_002E": "below_poverty",
    "B01003_001E": "total_population",
    "B03002_003E": "white_non_hispanic",
    "B03002_004E": "black_non_hispanic",
    "B03002_012E": "hispanic_latino",
}


def fetch_acs_data(api_key: str = None, years: list = None) -> pd.DataFrame | None:
    """
    Fetch ACS 1-year state-level estimates for a range of years.

    No API key is required — the Census API allows keyless access.
    Providing a key raises the rate limit but is optional.

    Note: ACS 1-year was not released in 2020 due to COVID-19 disruptions.
    2020 values are linearly interpolated from 2019 and 2021.

    Parameters
    ----------
    api_key : Census API key (optional; free at https://api.census.gov/data/key_signup.html)
    years   : list of years to fetch (default 2010–2022, skipping 2020)

    Returns
    -------
    DataFrame with demographic controls, or None on complete failure.
    """
    if years is None:
        years = [y for y in range(2010, 2023) if y != 2020]  # no ACS 1-yr in 2020

    frames = []
    for year in years:
        url = f"https://api.census.gov/data/{year}/acs/acs1"
        params = {"get": "NAME," + ",".join(_ACS_VARS.keys()), "for": "state:*"}
        if api_key:
            params["key"] = api_key
        try:
            r = requests.get(url, params=params, timeout=30)
            if r.status_code != 200 or r.text.strip().startswith("<"):
                print(f"  [WARN] {year}: API error (status {r.status_code})")
                continue

            data = r.json()
            df = pd.DataFrame(data[1:], columns=data[0])
            df["year"] = year
            df = df.rename(columns={"NAME": "state", "state": "state_fips",
                                    **_ACS_VARS})
            frames.append(df)
            print(f"  [OK] ACS {year}: {len(df)} states")

        except Exception as exc:
            print(f"  ⚠️  {year}: {exc}")

    if not frames:
        return None

    df = pd.concat(frames, ignore_index=True)

    # Calculate rates
    for col in list(_ACS_VARS.values()):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["poverty_rate"] = (df["below_poverty"] / df["poverty_universe"] * 100).round(2)
    df["pct_white"]    = (df["white_non_hispanic"] / df["total_population"] * 100).round(2)
    df["pct_black"]    = (df["black_non_hispanic"] / df["total_population"] * 100).round(2)
    df["pct_hispanic"] = (df["hispanic_latino"]    / df["total_population"] * 100).round(2)

    df["state_fips"] = df["state_fips"].str.zfill(2)
    df["median_household_income"] = df["median_household_income"].round(0)

    keep = ["state_fips", "year", "total_population", "median_household_income",
            "poverty_rate", "pct_white", "pct_black", "pct_hispanic"]
    df = df[keep].copy()

    # Interpolate 2020 (ACS 1-year not released due to COVID-19)
    all_fips = df["state_fips"].unique()
    rows_2020 = []
    for fips in all_fips:
        sub = df[df["state_fips"] == fips].set_index("year")
        if 2019 in sub.index and 2021 in sub.index:
            row_2020 = {"state_fips": fips, "year": 2020}
            for col in keep[2:]:
                row_2020[col] = round((sub.loc[2019, col] + sub.loc[2021, col]) / 2, 2)
            rows_2020.append(row_2020)
    if rows_2020:
        df = pd.concat([df, pd.DataFrame(rows_2020)], ignore_index=True)
        print(f"  [OK] ACS 2020: interpolated from 2019+2021 ({len(rows_2020)} states)")

    return df.sort_values(["state_fips", "year"]).reset_index(drop=True)


# ── Panel loader ──────────────────────────────────────────────────────────────

def load_panel(path: str = "../data/processed/analysis_panel.csv") -> pd.DataFrame:
    """Load the analysis panel and ensure state_fips is zero-padded."""
    df = pd.read_csv(path)
    df["state_fips"] = df["state_fips"].astype(str).str.zfill(2)
    return df
