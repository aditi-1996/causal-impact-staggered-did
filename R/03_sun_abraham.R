# ============================================================================
# 03 — Sun & Abraham (2021) Interaction-Weighted Estimator
# ============================================================================
#
# Implements the Sun & Abraham estimator which corrects TWFE bias by
# using interaction-weighted (IW) estimation.
#
# Key idea: Instead of a single treatment dummy, interact treatment
# with cohort indicators and use never-treated/last-treated as controls.
#
# Reference: Sun & Abraham (2021), Journal of Econometrics
# Implementation via fixest::sunab()
# ============================================================================

library(fixest)
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

panel <- panel %>%
  mutate(
    state_id = as.numeric(as.factor(state_fips)),
    # For sunab: cohort_group must be numeric
    # Never-treated states need a large value (Inf or very large year)
    cohort_sunab = ifelse(cohort_group == 0, Inf, cohort_group)
  )

cat("Panel loaded:", nrow(panel), "rows\n")

# ── Identify outcomes ────────────────────────────────────────────────────────
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

cat("Outcomes:", paste(outcome_vars, collapse = ", "), "\n")

# ── Controls ─────────────────────────────────────────────────────────────────
control_candidates <- c(
  "median_household_income", "poverty_rate", "total_population",
  "pct_white", "pct_black", "pct_hispanic"
)
control_vars <- control_candidates[
  control_candidates %in% names(panel) &
  sapply(control_candidates, function(v) sum(!is.na(panel[[v]])) > 200)
]

# ── Sun & Abraham Estimation ─────────────────────────────────────────────────
sa_results <- list()

for (outcome in outcome_vars) {
  cat("\n", strrep("=", 60), "\n")
  cat("Sun & Abraham:", outcome, "\n")
  cat(strrep("=", 60), "\n")
  
  df <- panel %>% filter(!is.na(.data[[outcome]]))
  
  # Build formula: sunab(cohort, time) replaces treatment variable
  ctrl <- setdiff(control_vars, outcome)
  ctrl_complete <- ctrl[sapply(ctrl, function(v) sum(is.na(df[[v]])) < nrow(df) * 0.1)]
  ctrl_str <- if (length(ctrl_complete) > 0) paste("+", paste(ctrl_complete, collapse = " + ")) else ""
  
  fml <- as.formula(paste0(
    outcome, " ~ sunab(cohort_sunab, year) ", ctrl_str, " | state_id + year"
  ))
  
  tryCatch({
    est <- feols(fml, data = df, cluster = ~state_id)
    sa_results[[outcome]] <- est
    
    cat("\nSummary:\n")
    print(summary(est))
    
    # ── Aggregate ATT ──
    agg <- summary(est, agg = "ATT")
    cat("\nAggregated ATT:\n")
    print(agg)
    
    # ── Event Study Plot ──
    coef_names <- names(coef(est))
    # Extract event time from coefficient names
    et_pattern <- grepl("year::", coef_names)
    
    if (sum(et_pattern) > 0) {
      coef_df <- data.frame(
        event_time = as.numeric(gsub(".*::", "", coef_names[et_pattern])),
        estimate = coef(est)[et_pattern],
        ci_low = confint(est)[et_pattern, 1],
        ci_high = confint(est)[et_pattern, 2]
      )
      
      # Add reference period
      coef_df <- rbind(coef_df, data.frame(
        event_time = -1, estimate = 0, ci_low = 0, ci_high = 0
      ))
      coef_df <- coef_df[order(coef_df$event_time), ]
      
      p <- ggplot(coef_df, aes(x = event_time, y = estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.5) +
        geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "#238b45", alpha = 0.2) +
        geom_point(color = "#238b45", size = 2.5) +
        geom_line(color = "#238b45", linewidth = 0.8) +
        labs(
          title = paste("Sun & Abraham Event Study:", gsub("_", " ", outcome)),
          subtitle = "Interaction-weighted estimator | State + Year FE | Clustered SEs",
          x = "Years Relative to Medicaid Expansion",
          y = "Estimated Effect (95% CI)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        )
      
      ggsave(
        file.path(OUTPUT_DIR, paste0("sa_event_study_", outcome, ".png")),
        plot = p, width = 10, height = 6, dpi = 150
      )
      cat("Saved: figures/sa_event_study_", outcome, ".png\n")
    }
    
  }, error = function(e) {
    cat("⚠️  Error:", conditionMessage(e), "\n")
  })
}

# ── Compare TWFE vs Sun-Abraham ──────────────────────────────────────────────
cat("\n\n", strrep("=", 60), "\n")
cat("TWFE vs Sun & Abraham Comparison\n")
cat(strrep("=", 60), "\n")

for (outcome in outcome_vars) {
  if (outcome %in% names(sa_results)) {
    df <- panel %>% filter(!is.na(.data[[outcome]]))
    
    ctrl <- setdiff(control_vars, outcome)
    ctrl_complete <- ctrl[sapply(ctrl, function(v) sum(is.na(df[[v]])) < nrow(df) * 0.1)]
    ctrl_str <- if (length(ctrl_complete) > 0) paste("+", paste(ctrl_complete, collapse = " + ")) else ""
    
    # Standard TWFE
    fml_twfe <- as.formula(paste0(
      outcome, " ~ post_expansion ", ctrl_str, " | state_id + year"
    ))
    
    tryCatch({
      twfe <- feols(fml_twfe, data = df, cluster = ~state_id)
      
      cat("\n", outcome, ":\n")
      cat("  TWFE estimate:          ", round(coef(twfe)["post_expansion"], 4), "\n")
      cat("  Sun-Abraham aggregate:  see summary above\n")
      cat("  Difference suggests:    ",
          "potential bias in TWFE due to treatment effect heterogeneity\n")
    }, error = function(e) {
      cat("⚠️  TWFE comparison error:", conditionMessage(e), "\n")
    })
  }
}

cat("\n✅ Sun & Abraham analysis complete.\n")
