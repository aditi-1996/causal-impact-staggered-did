# ============================================================================
# 04 — TWFE Diagnostics: Goodman-Bacon Decomposition
# ============================================================================
#
# The Goodman-Bacon (2021) decomposition shows how the standard TWFE
# estimator is a weighted average of all possible 2x2 DiD comparisons,
# including potentially biased "bad comparisons" that use already-treated
# units as controls.
#
# This script:
#   1. Runs the Bacon decomposition to identify sub-estimates
#   2. Visualizes the decomposition
#   3. Flags potential bias from timing-group comparisons
#
# Reference: Goodman-Bacon (2021), Journal of Econometrics
# ============================================================================

library(bacondecomp)
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
  mutate(state_id = as.numeric(as.factor(state_fips)))

cat("Panel loaded:", nrow(panel), "rows\n")

# ── Identify outcomes ────────────────────────────────────────────────────────
outcome_candidates <- c(
  "allcause_age_adj_rate", "allcause_crude_rate",
  "diabetes_age_adj_rate", "diabetes_crude_rate",
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

# ── Bacon Decomposition ─────────────────────────────────────────────────────
bacon_results <- list()

for (outcome in outcome_vars) {
  cat("\n", strrep("=", 60), "\n")
  cat("Bacon Decomposition:", outcome, "\n")
  cat(strrep("=", 60), "\n")
  
  # bacondecomp requires a balanced panel with no NAs
  df <- panel %>%
    filter(!is.na(.data[[outcome]])) %>%
    select(state_id, year, all_of(outcome), post_expansion) %>%
    # Ensure balanced panel
    group_by(state_id) %>%
    filter(n() == n_distinct(panel$year)) %>%
    ungroup()
  
  if (nrow(df) == 0) {
    cat("⚠️  No balanced panel available for", outcome, "\n")
    next
  }
  
  tryCatch({
    # Run decomposition
    bacon_out <- bacon(
      as.formula(paste(outcome, "~ post_expansion")),
      data = as.data.frame(df),
      id_var = "state_id",
      time_var = "year"
    )
    
    bacon_results[[outcome]] <- bacon_out
    
    # ── Summary ──
    cat("\nDecomposition summary:\n")
    type_summary <- bacon_out %>%
      group_by(type) %>%
      summarise(
        n_comparisons = n(),
        avg_estimate = round(mean(estimate), 4),
        total_weight = round(sum(weight), 4),
        weighted_estimate = round(sum(estimate * weight) / sum(weight), 4),
        .groups = "drop"
      )
    print(as.data.frame(type_summary))
    
    # Overall TWFE estimate (weighted average)
    twfe_estimate <- sum(bacon_out$estimate * bacon_out$weight)
    cat("\nOverall TWFE estimate (weighted avg):", round(twfe_estimate, 4), "\n")
    
    # ── Bacon Plot ──
    p <- ggplot(bacon_out, aes(x = weight, y = estimate, color = type, shape = type)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_hline(yintercept = twfe_estimate, linetype = "dashed", color = "black", linewidth = 0.5) +
      annotate("text", x = max(bacon_out$weight) * 0.8, y = twfe_estimate,
               label = paste("TWFE =", round(twfe_estimate, 3)),
               vjust = -1, size = 3.5) +
      geom_hline(yintercept = 0, color = "gray50", linewidth = 0.3) +
      scale_color_manual(values = c(
        "Earlier vs Later Treated" = "#e6550d",
        "Later vs Earlier Treated" = "#d94701",
        "Treated vs Untreated" = "#2171b5"
      )) +
      labs(
        title = paste("Goodman-Bacon Decomposition:", gsub("_", " ", outcome)),
        subtitle = "Each point is a 2x2 DiD comparison; size ∝ weight in TWFE estimate",
        x = "Weight in TWFE Estimate",
        y = "2x2 DiD Estimate",
        color = "Comparison Type",
        shape = "Comparison Type"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      file.path(OUTPUT_DIR, paste0("bacon_decomp_", outcome, ".png")),
      plot = p, width = 10, height = 7, dpi = 150
    )
    cat("Saved: figures/bacon_decomp_", outcome, ".png\n")
    
    # ── Diagnose bias ──
    bad_comparisons <- bacon_out %>%
      filter(type %in% c("Earlier vs Later Treated", "Later vs Earlier Treated"))
    
    if (nrow(bad_comparisons) > 0) {
      bad_weight <- sum(bad_comparisons$weight)
      cat("\n⚠️  'Bad comparison' weight:", round(bad_weight * 100, 1), "% of TWFE estimate\n")
      cat("This suggests", if (bad_weight > 0.3) "SUBSTANTIAL" else "moderate",
          "potential for bias.\n")
      cat("→ Use Callaway-Sant'Anna or Sun-Abraham estimators for robust estimates.\n")
    } else {
      cat("\n✅ No timing-based comparisons found — TWFE may be unbiased.\n")
    }
    
  }, error = function(e) {
    cat("⚠️  Error:", conditionMessage(e), "\n")
    cat("bacondecomp requires a balanced panel. Check for missing state-year observations.\n")
  })
}

# ── Save decomposition tables ────────────────────────────────────────────────
if (length(bacon_results) > 0) {
  for (outcome in names(bacon_results)) {
    write_csv(
      bacon_results[[outcome]],
      file.path(TABLE_DIR, paste0("bacon_decomp_", outcome, ".csv"))
    )
  }
  cat("\nSaved decomposition tables to tables/\n")
}

cat("\n✅ Goodman-Bacon diagnostics complete.\n")
