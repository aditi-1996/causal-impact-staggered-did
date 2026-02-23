# ============================================================================
# 02 — Event Study with fixest
# ============================================================================
#
# Fast event study estimation using the fixest package.
# Produces clean event study plots with confidence intervals.
#
# fixest advantages:
#   - Extremely fast fixed effects estimation
#   - Built-in event study support via i() function
#   - Easy clustered standard errors
#   - Publication-quality coefficient plots
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

# Create numeric state ID for fixed effects
panel <- panel %>%
  mutate(
    state_id = as.numeric(as.factor(state_fips)),
    # Bin event time at endpoints for clean estimation
    event_time_binned = case_when(
      is.na(event_time) ~ NA_real_,  # Never-treated
      event_time < -4 ~ -4,
      event_time > 8 ~ 8,
      TRUE ~ event_time
    )
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

# ── Control formula ──────────────────────────────────────────────────────────
control_candidates <- c(
  "median_household_income", "poverty_rate", "total_population",
  "pct_white", "pct_black", "pct_hispanic"
)

control_vars <- control_candidates[
  control_candidates %in% names(panel) &
  sapply(control_candidates, function(v) sum(!is.na(panel[[v]])) > 200)
]

# ── TWFE Event Study with fixest ─────────────────────────────────────────────
es_results <- list()

for (outcome in outcome_vars) {
  cat("\n", strrep("=", 60), "\n")
  cat("Event Study (fixest):", outcome, "\n")
  cat(strrep("=", 60), "\n")
  
  df <- panel %>% filter(!is.na(.data[[outcome]]))
  
  # Build formula with event study interaction
  # i(event_time_binned, ref = -1) creates dummies with t=-1 as reference
  ctrl_str <- ""
  ctrl <- setdiff(control_vars, outcome)
  ctrl_complete <- ctrl[sapply(ctrl, function(v) sum(is.na(df[[v]])) < nrow(df) * 0.1)]
  if (length(ctrl_complete) > 0) {
    ctrl_str <- paste("+", paste(ctrl_complete, collapse = " + "))
  }
  
  fml <- as.formula(paste0(
    outcome, " ~ i(event_time_binned, ref = -1) ", ctrl_str,
    " | state_id + year"
  ))
  
  tryCatch({
    # Estimate
    est <- feols(fml, data = df, cluster = ~state_id)
    es_results[[outcome]] <- est
    
    cat("\nSummary:\n")
    print(summary(est))
    
    # ── Event Study Plot ──
    p <- iplot(est, main = paste("Event Study:", outcome),
               xlab = "Years Relative to Expansion",
               ylab = "Estimated Effect")
    
    # Save using ggplot for better control
    coef_df <- data.frame(
      event_time = as.numeric(gsub("event_time_binned::", "", names(coef(est)))),
      estimate = coef(est),
      ci_low = confint(est)[, 1],
      ci_high = confint(est)[, 2]
    )
    # Filter to event time coefficients only
    coef_df <- coef_df[!is.na(coef_df$event_time), ]
    
    # Add reference period
    coef_df <- rbind(coef_df, data.frame(
      event_time = -1, estimate = 0, ci_low = 0, ci_high = 0
    ))
    coef_df <- coef_df[order(coef_df$event_time), ]
    
    p_gg <- ggplot(coef_df, aes(x = event_time, y = estimate)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.5) +
      geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "#2171b5", alpha = 0.2) +
      geom_point(color = "#2171b5", size = 2.5) +
      geom_line(color = "#2171b5", linewidth = 0.8) +
      annotate("text", x = -2.5, y = min(coef_df$ci_low) * 0.9,
               label = "Pre-treatment", color = "gray50", fontface = "italic") +
      annotate("text", x = 4, y = min(coef_df$ci_low) * 0.9,
               label = "Post-treatment", color = "gray50", fontface = "italic") +
      labs(
        title = paste("Event Study (TWFE):", gsub("_", " ", outcome)),
        subtitle = "State + Year FE | Clustered SEs",
        x = "Years Relative to Medicaid Expansion",
        y = "Estimated Effect (95% CI)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      file.path(OUTPUT_DIR, paste0("fixest_event_study_", outcome, ".png")),
      plot = p_gg, width = 10, height = 6, dpi = 150
    )
    cat("Saved: figures/fixest_event_study_", outcome, ".png\n")
    
  }, error = function(e) {
    cat("⚠️  Error:", conditionMessage(e), "\n")
  })
}

# ── Comparison table ─────────────────────────────────────────────────────────
if (length(es_results) > 1) {
  cat("\n\n", strrep("=", 60), "\n")
  cat("Side-by-side comparison\n")
  cat(strrep("=", 60), "\n")
  etable(es_results, se = "cluster")
}

cat("\n✅ fixest event study analysis complete.\n")
