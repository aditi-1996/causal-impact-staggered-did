# ============================================================================
# 01 — Callaway & Sant'Anna (2021) Difference-in-Differences
# ============================================================================
#
# This script implements the modern DiD estimator that avoids the bias
# in standard TWFE with staggered treatment adoption.
#
# Key advantages over TWFE:
#   - Never uses already-treated units as controls
#   - Estimates group-time specific ATTs
#   - Aggregates to overall, dynamic, and group-specific effects
#
# Reference: Callaway & Sant'Anna (2021), Journal of Econometrics
# ============================================================================

library(did)
library(dplyr)
library(ggplot2)
library(readr)

# ── Settings ──────────────────────────────────────────────────────────────────
OUTPUT_DIR <- "../figures"
TABLE_DIR <- "../tables"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load Data ─────────────────────────────────────────────────────────────────
panel <- read_csv("../data/processed/analysis_panel.csv", show_col_types = FALSE)

cat("Panel loaded:", nrow(panel), "rows,", n_distinct(panel$state), "states\n")
cat("Years:", min(panel$year), "–", max(panel$year), "\n")

# ── Prepare for did package ──────────────────────────────────────────────────
# The did package requires:
#   - yname: outcome variable
#   - tname: time variable
#   - idname: unit identifier (numeric)
#   - gname: group variable (first treatment period, 0 for never-treated)

# Create numeric state ID; coerce time/group vars to integer for did package
panel <- panel %>%
  mutate(
    state_id     = as.integer(as.factor(state_fips)),
    year         = as.integer(year),
    cohort_group = as.integer(cohort_group)
  )

# Verify cohort_group coding (0 = never treated)
cat("\nCohort groups:\n")
print(table(panel$cohort_group))

# ── Identify available outcome variables ─────────────────────────────────────
outcome_candidates <- c(
  "allcause_age_adj_rate", "allcause_crude_rate",
  "diabetes_age_adj_rate", "diabetes_crude_rate",
  "maternal_age_adj_rate", "maternal_crude_rate",
  "diabetes_prevalence", "low_birth_weight_pct",
  "poverty_rate", "median_household_income"
)

outcome_vars <- outcome_candidates[
  outcome_candidates %in% names(panel) &
  sapply(outcome_candidates, function(v) {
    if (v %in% names(panel)) sum(!is.na(panel[[v]])) > 100 else FALSE
  })
]

cat("\nAvailable outcomes:", paste(outcome_vars, collapse = ", "), "\n")

if (length(outcome_vars) == 0) {
  stop("No outcome variables available. Download CDC data and re-run cleaning notebook.")
}

# ── Control variables ────────────────────────────────────────────────────────
control_candidates <- c(
  "median_household_income", "poverty_rate", "total_population",
  "pct_white", "pct_black", "pct_hispanic"
)

control_vars <- control_candidates[
  control_candidates %in% names(panel) &
  sapply(control_candidates, function(v) {
    if (v %in% names(panel)) sum(!is.na(panel[[v]])) > 200 else FALSE
  })
]

# Remove outcome from controls if it appears
use_controls <- length(control_vars) > 0

cat("Control variables:", if (use_controls) paste(control_vars, collapse = ", ") else "NONE", "\n")

# ── Run Callaway-Sant'Anna for each outcome ──────────────────────────────────
cs_results <- list()

for (outcome in outcome_vars) {
  cat("\n", strrep("=", 60), "\n")
  cat("Callaway-Sant'Anna:", outcome, "\n")
  cat(strrep("=", 60), "\n")
  
  # Keep only rows with non-missing outcome; minimal columns for att_gt
  # Note: covariates cause aggte to crash in did 2.1.2 — run unconditional.
  # Controlled estimates are produced by run_cohort_did() in the Python pipeline.
  keep_cols <- c("state_id", "year", "cohort_group", outcome)
  df <- as.data.frame(panel[!is.na(panel[[outcome]]), keep_cols])
  xformla <- NULL

  tryCatch({
    # Maternal outcomes have too few never-treated obs after NA filtering
    ctrl_grp <- if (grepl("maternal", outcome)) "notyettreated" else "nevertreated"

    # Estimate group-time ATTs
    cs_out <- att_gt(
      yname = outcome,
      tname = "year",
      idname = "state_id",
      gname = "cohort_group",
      xformla = xformla,
      data = df,
      control_group = ctrl_grp,
      anticipation = 0,
      est_method = "reg",
      base_period = "varying",
      clustervars = "state_id",
      print_details = FALSE
    )
    
    cs_results[[outcome]] <- cs_out
    
    # ── Overall ATT ──
    agg_overall <- aggte(cs_out, type = "simple", na.rm = TRUE)
    cat("\nOverall ATT:", round(agg_overall$overall.att, 4), "\n")
    cat("SE:", round(agg_overall$overall.se, 4), "\n")
    cat("95% CI: [", round(agg_overall$overall.att - 1.96 * agg_overall$overall.se, 4),
        ",", round(agg_overall$overall.att + 1.96 * agg_overall$overall.se, 4), "]\n")

    # ── Dynamic (Event Study) Aggregation ──
    agg_dynamic <- aggte(cs_out, type = "dynamic", na.rm = TRUE)

    # Plot event study
    p <- ggdid(agg_dynamic) +
      theme_minimal(base_size = 12) +
      labs(
        title = paste("Callaway-Sant'Anna Event Study:", outcome),
        subtitle = "Group-time ATTs aggregated by event time",
        x = "Years Relative to Expansion",
        y = "Estimated ATT"
      ) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
    
    ggsave(
      file.path(OUTPUT_DIR, paste0("cs_event_study_", outcome, ".png")),
      plot = p, width = 10, height = 6, dpi = 150
    )
    cat("Saved: figures/cs_event_study_", outcome, ".png\n")
    
    # ── Group-specific ATTs ──
    agg_group <- aggte(cs_out, type = "group", na.rm = TRUE)
    cat("\nGroup-specific ATTs:\n")
    group_df <- data.frame(
      group = agg_group$egt,
      att = round(agg_group$att.egt, 4),
      se = round(agg_group$se.egt, 4)
    )
    print(group_df)
    
  }, error = function(e) {
    cat("\n⚠️  Error:", conditionMessage(e), "\n")
  })
}

# ── Save summary ─────────────────────────────────────────────────────────────
if (length(cs_results) > 0) {
  summary_rows <- lapply(names(cs_results), function(outcome) {
    agg <- aggte(cs_results[[outcome]], type = "simple", na.rm = TRUE)
    data.frame(
      outcome = outcome,
      method = "Callaway-Sant'Anna",
      att = round(agg$overall.att, 4),
      se = round(agg$overall.se, 4),
      ci_lower = round(agg$overall.att - 1.96 * agg$overall.se, 4),
      ci_upper = round(agg$overall.att + 1.96 * agg$overall.se, 4)
    )
  })
  
  summary_df <- do.call(rbind, summary_rows)
  write_csv(summary_df, file.path(TABLE_DIR, "cs_results_summary.csv"))
  cat("\n\nSaved: tables/cs_results_summary.csv\n")
  print(summary_df)
}

cat("\n✅ Callaway-Sant'Anna analysis complete.\n")
