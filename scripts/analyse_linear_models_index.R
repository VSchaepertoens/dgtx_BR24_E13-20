#===============================================================================
# GLYCOSYLATION INDEX ANALYSIS - SPLINE MODEL FITTING AND INTERPRETATION
#===============================================================================

# Load required libraries
library(tidyverse)
library(splines)
library(here)
library(multcomp)
library(gridExtra)
#source("scripts/plotting_utils.R", local = FALSE)

# This works both interactively and when knitting
project_root <- here::here()

# Source plotting utils
source(file.path(project_root, "scripts", "plotting_utils.R"), local = FALSE)

#===============================================================================
# CONFIGURATION - CHANGE THIS LINE FOR DIFFERENT INDICES
#===============================================================================

# Options: "galactosylation", "fucosylation", "glycation"
index_type <- "fucosylation"  # <-- CHANGE THIS LINE FOR DIFFERENT INDICES

#===============================================================================
# 1. DATA LOADING AND PREPARATION
#===============================================================================

# Set file path and display names based on index_type
# input_file_path <- here::here("analysis", paste0("charrun_E13-E20_", index_type, "index_V01_20260114_VS.RData"))
input_file_path <- file.path(project_root, "analysis", 
                             paste0("charrun_E13-E20_", index_type, "index_V01_20260114_VS.RData"))
load(file = input_file_path)

# Determine display name and value column
display_name <- paste0(toupper(substr(index_type, 1, 1)), substring(index_type, 2))
value_col <- ifelse(index_type == "galactosylation", "GI",
                    ifelse(index_type == "fucosylation", "GI",
                           ifelse(index_type == "glycation", "GI")))

# Set color scheme for conditions
color_mapping_condition <- c(
  "CT" = "#E6641E",
  "TS" = "#4B288C"
)

# Set reference level for condition (CT as baseline)
gi_summary$condition <- factor(gi_summary$condition)
gi_summary$condition <- relevel(gi_summary$condition, ref = "CT")
gi_summary$experiment <- factor(gi_summary$experiment)

# Create formula dynamically
formula_spline <- as.formula(paste0(value_col, " ~ ns(tp, df = 3) * condition"))
formula_lin <- as.formula(paste0(value_col, " ~ tp * condition"))
formula_quad <- as.formula(paste0(value_col, " ~ (tp + I(tp^2)) * condition"))

output_fig_dir <- here::here("figures", paste0(index_type, "_index"))
if (!dir.exists(output_fig_dir)) dir.create(output_fig_dir, recursive = TRUE)

#===============================================================================
# 2. EXPLORATORY DATA ANALYSIS
#===============================================================================

cat("\n========== CHECKING LINEARITY ASSUMPTION ==========\n")

# Visual inspection
p_linearity <- ggplot(gi_summary, aes(tp, .data[[value_col]], color = condition)) +
  geom_point(aes(shape = experiment), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  scale_color_manual(values = color_mapping_condition, name = "Condition") +
  scale_shape_manual(values = 1:nlevels(gi_summary$experiment), name = "Experiment") +
  labs(title = paste(display_name, "- Linearity Assessment"),
       x = "Time [h]", y = paste0(display_name, " Index [%]")) +
  theme_minimal()

print(p_linearity)

# Formal tests for linearity
fit_lin <- lm(formula_lin, data = gi_summary)
fit_quad <- lm(formula_quad, data = gi_summary)
fit_spline <- lm(formula_spline, data = gi_summary)

cat("\nModel Comparison:\n")
cat("- Linear vs Quadratic: ")
print(anova(fit_lin, fit_quad))
cat("- Linear vs Spline (df=3): ")
print(anova(fit_lin, fit_spline))

# 2.2 Check for outliers
p_outliers <- ggplot(gi_summary, aes(y = .data[[value_col]], x = condition, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = color_mapping_condition, guide = "none") +
  labs(title = paste(display_name, "- Outlier Check"), 
       x = "Condition", y = paste0(display_name, " Index [%]")) +
  theme_minimal()

print(p_outliers)

# Save outlier plot
# ggsave(file.path(output_fig_dir, paste0(index_type, "_outlier_check.png")), 
#        p_outliers, width = 6, height = 5, dpi = 300, bg = "white")
# ggsave(file.path(output_fig_dir, paste0(index_type, "_outlier_check.pdf")), 
#        p_outliers, width = 6, height = 5)

#===============================================================================
# 3. SPLINE MODEL FITTING
#===============================================================================

cat("\n========== SPLINE MODEL FIT ==========\n")
print(summary(fit_spline))

#===============================================================================
# 4. STATISTICAL INFERENCE USING CONTRAST MATRICES
#===============================================================================

# Extract coefficient names
coef_names <- names(coef(fit_spline))
spline_terms <- grep("^ns\\(tp, df = 3\\)[0-9]$", coef_names, value = TRUE)
interaction_terms <- grep(":conditionTS", coef_names, value = TRUE)

#---------------------------------------------------------------------------
# 4.1 Time Effect in CT (Reference Condition)
#---------------------------------------------------------------------------
L_time_CT <- matrix(0, nrow = length(spline_terms), ncol = length(coef_names))
colnames(L_time_CT) <- coef_names
for (i in seq_along(spline_terms)) {
  L_time_CT[i, spline_terms[i]] <- 1
}
time_CT_test <- glht(fit_spline, linfct = L_time_CT)
time_CT_result <- summary(time_CT_test, test = Ftest())

#---------------------------------------------------------------------------
# 4.2 Time Effect in TS Condition
#---------------------------------------------------------------------------
L_time_TS <- matrix(0, nrow = length(spline_terms), ncol = length(coef_names))
colnames(L_time_TS) <- coef_names
for (i in seq_along(spline_terms)) {
  L_time_TS[i, spline_terms[i]] <- 1
  L_time_TS[i, interaction_terms[i]] <- 1
}
time_TS_test <- glht(fit_spline, linfct = L_time_TS)
time_TS_result <- summary(time_TS_test, test = Ftest())

#---------------------------------------------------------------------------
# 4.3 Condition Effect (Average Difference)
#---------------------------------------------------------------------------
L_condition <- matrix(0, nrow = 1, ncol = length(coef_names))
colnames(L_condition) <- coef_names
L_condition[1, "conditionTS"] <- 1
condition_test <- glht(fit_spline, linfct = L_condition)
condition_result <- summary(condition_test, test = Ftest())

#---------------------------------------------------------------------------
# 4.4 Interaction Effect (Time × Condition)
#---------------------------------------------------------------------------
L_interaction <- matrix(0, nrow = length(interaction_terms), ncol = length(coef_names))
colnames(L_interaction) <- coef_names
for (i in seq_along(interaction_terms)) {
  L_interaction[i, interaction_terms[i]] <- 1
}
interaction_test <- glht(fit_spline, linfct = L_interaction)
interaction_result <- summary(interaction_test, test = Ftest())

#===============================================================================
# 5. EFFECT SIZE METRICS (cT, CV, cDT)
#===============================================================================

# Function to calculate spline-based metrics
calculate_spline_metrics <- function(fit, data, condition_val = NULL) {
  tp_range <- range(data$tp)
  tp_grid <- seq(tp_range[1], tp_range[2], length.out = 1000)
  
  if (is.null(condition_val)) {
    # Average across both conditions
    pred_CT <- predict(fit, data.frame(tp = tp_grid, condition = "CT"))
    pred_TS <- predict(fit, data.frame(tp = tp_grid, condition = "TS"))
    pred <- (pred_CT + pred_TS) / 2
  } else {
    pred <- predict(fit, data.frame(tp = tp_grid, condition = condition_val))
  }
  
  list(
    cv = sd(pred) / mean(pred) * 100,
    cumulative_travel = sum(abs(diff(pred))),
    range_change = max(pred) - min(pred),
    predictions = data.frame(tp = tp_grid, pred = pred)
  )
}

# Calculate metrics for each component
metrics_CT <- calculate_spline_metrics(fit_spline, gi_summary, "CT")
metrics_TS <- calculate_spline_metrics(fit_spline, gi_summary, "TS")
metrics_avg <- calculate_spline_metrics(fit_spline, gi_summary, NULL)

# Calculate interaction-specific metrics
tp_grid <- seq(min(gi_summary$tp), max(gi_summary$tp), length.out = 1000)
pred_CT <- predict(fit_spline, data.frame(tp = tp_grid, condition = "CT"))
pred_TS <- predict(fit_spline, data.frame(tp = tp_grid, condition = "TS"))
diff_trajectory <- pred_TS - pred_CT

cDT <- sum(abs(diff(diff_trajectory)))        # Cumulative Differential Travel
avg_abs_diff <- mean(abs(diff_trajectory))    # Average absolute difference
max_divergence <- max(abs(diff_trajectory))   # Maximum divergence

#===============================================================================
# 6. COMPREHENSIVE RESULTS SUMMARY
#===============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║           ", toupper(index_type), "INDEX ANALYSIS RESULTS                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")

# 6.1 Statistical significance table
results_table <- data.frame(
  Effect = c("Time (CT condition)", "Time (TS condition)", "Condition", "Time × Condition"),
  DF = c(length(spline_terms), length(spline_terms), 1, length(interaction_terms)),
  F_statistic = c(
    time_CT_result$test$fstat[1],
    time_TS_result$test$fstat[1],
    condition_result$test$fstat[1],
    interaction_result$test$fstat[1]
  ),
  P_value = c(
    time_CT_result$test$pvalue[1],
    time_TS_result$test$pvalue[1],
    condition_result$test$pvalue[1],
    interaction_result$test$pvalue[1]
  )
)

cat("\n┌───────────────────────────── STATISTICAL SIGNIFICANCE ──────────────────────────┐\n")
print(results_table, row.names = FALSE)
cat("└─────────────────────────────────────────────────────────────────────────────────┘\n")

# 6.2 Effect size metrics table
metrics_table <- data.frame(
  Metric = c(paste0("Cumulative Travel (cT) [", value_col, " units]"), 
             "Coefficient of Variation (CV) [%]",
             paste0("Range Change [", value_col, " units]")),
  CT = c(round(metrics_CT$cumulative_travel, 3),
         round(metrics_CT$cv, 1),
         round(metrics_CT$range_change, 3)),
  TS = c(round(metrics_TS$cumulative_travel, 3),
         round(metrics_TS$cv, 1),
         round(metrics_TS$range_change, 3)),
  Ratio_TS_CT = c(round(metrics_TS$cumulative_travel / metrics_CT$cumulative_travel, 2),
                  round(metrics_TS$cv / metrics_CT$cv, 2),
                  round(metrics_TS$range_change / metrics_CT$range_change, 2))
)

cat("\n┌───────────────────────────── EFFECT SIZE METRICS ───────────────────────────────┐\n")
print(metrics_table, row.names = FALSE)
cat("└─────────────────────────────────────────────────────────────────────────────────┘\n")

interaction_metrics <- data.frame(
  Metric = c("Cumulative Differential Travel (cDT)", 
             "Average Absolute Difference",
             "Maximum Divergence"),
  Value = c(round(cDT, 3), round(avg_abs_diff, 3), round(max_divergence, 3)),
  Units = rep(paste0(value_col, " units"), 3)
)

cat("\n┌───────────────────────────── INTERACTION METRICS ───────────────────────────────┐\n")
print(interaction_metrics, row.names = FALSE)
cat("└─────────────────────────────────────────────────────────────────────────────────┘\n")

#===============================================================================
# 7. INTERPRETATION
#===============================================================================

cat("\n╔══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                                 INTERPRETATION                                ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")

# Helper function for formatting p-values
format_p <- function(p) {
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

# Time effect interpretation
cat("\n▶ TIME EFFECT:\n")
if (time_CT_result$test$pvalue[1] < 0.05) {
  cat("  • CT condition: Significant temporal change (F =", 
      round(time_CT_result$test$fstat[1], 2), ", p", format_p(time_CT_result$test$pvalue[1]), ")\n")
} else {
  cat("  • CT condition: No significant temporal change (p", format_p(time_CT_result$test$pvalue[1]), ")\n")
}

if (time_TS_result$test$pvalue[1] < 0.05) {
  cat("  • TS condition: Significant temporal change (F =", 
      round(time_TS_result$test$fstat[1], 2), ", p", format_p(time_TS_result$test$pvalue[1]), ")\n")
} else {
  cat("  • TS condition: No significant temporal change (p", format_p(time_TS_result$test$pvalue[1]), ")\n")
}

cat("  • Effect size: CT cumulative travel =", round(metrics_CT$cumulative_travel, 3), 
    paste0(value_col, " units, TS ="), round(metrics_TS$cumulative_travel, 3), paste0(value_col, " units\n"))

# Condition effect interpretation
cat("\n▶ CONDITION EFFECT (Average difference across all timepoints):\n")
if (condition_result$test$pvalue[1] < 0.05) {
  cat("  • Significant difference between CT and TS (p", format_p(condition_result$test$pvalue[1]), ")\n")
} else {
  cat("  • No significant difference between CT and TS (p", format_p(condition_result$test$pvalue[1]), ")\n")
}
cat("  • Average absolute difference:", round(avg_abs_diff, 3), paste0(value_col, " units\n"))

# Interaction effect interpretation
cat("\n▶ INTERACTION EFFECT (Difference in temporal patterns):\n")
if (interaction_result$test$pvalue[1] < 0.05) {
  cat("  • Significant interaction: temporal patterns differ between conditions\n")
  cat("    (F =", round(interaction_result$test$fstat[1], 2), ", p", 
      format_p(interaction_result$test$pvalue[1]), ")\n")
} else {
  cat("  • No significant interaction: temporal patterns are similar between conditions\n")
  cat("    (F =", round(interaction_result$test$fstat[1], 2), ", p", 
      format_p(interaction_result$test$pvalue[1]), ")\n")
}
cat("  • Cumulative differential travel (cDT):", round(cDT, 3), paste0(value_col, " units\n"))

# Overall conclusion
cat("\n▶ OVERALL CONCLUSION:\n")
if (time_CT_result$test$pvalue[1] < 0.05 && time_TS_result$test$pvalue[1] < 0.05) {
  if (interaction_result$test$pvalue[1] > 0.05) {
    cat("  Both conditions show significant temporal changes with similar trajectories.\n")
    cat("  The ", tolower(display_name), " index changes substantially over time (cT ≈", 
        round((metrics_CT$cumulative_travel + metrics_TS$cumulative_travel)/2, 3), 
        paste0(value_col, " units average travel).\n"))
  } else {
    cat("  Both conditions show significant temporal changes, but with different trajectories.\n")
    cat("  The differential travel (cDT =", round(cDT, 3), 
        paste0(value_col, " units) quantifies the divergence between conditions.\n"))
  }
} else if (time_CT_result$test$pvalue[1] < 0.05) {
  cat("  Only the CT condition shows significant temporal changes.\n")
} else if (time_TS_result$test$pvalue[1] < 0.05) {
  cat("  Only the TS condition shows significant temporal changes.\n")
} else {
  cat("  Neither condition shows significant temporal changes.\n")
}

cat("\n")

#===============================================================================
# 7.5 SAVE RESULTS TO CSV FILES
#===============================================================================

# Create output directory if it doesn't exist
output_dir <- here::here("analysis", paste0(index_type, "_index"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Base filename with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
base_filename <- file.path(output_dir, paste0(index_type, "_spline_analysis_", timestamp))

#---------------------------------------------------------------------------
# Table 1: Combined Summary (All Results in One Table)
#---------------------------------------------------------------------------
combined_summary <- data.frame(
  Category = c(
    "STATISTICS", "", "", "",
    "CONDITION METRICS", "", "", "", "", "",
    "INTERACTION METRICS", "", "", ""
  ),
  Parameter = c(
    "Time effect (CT)", "Time effect (TS)", "Condition effect", "Interaction effect",
    "CT cumulative travel (cT)", "TS cumulative travel (cT)", "Ratio TS/CT",
    "CT CV", "TS CV", "CT range change", "TS range change",
    "Cumulative diff travel (cDT)", "Avg absolute difference", "Max divergence"
  ),
  Value = c(
    sprintf("F=%.2f, p=%s", time_CT_result$test$fstat[1], format(time_CT_result$test$pvalue[1], scientific = TRUE, digits = 3)),
    sprintf("F=%.2f, p=%s", time_TS_result$test$fstat[1], format(time_TS_result$test$pvalue[1], scientific = TRUE, digits = 3)),
    sprintf("F=%.2f, p=%.3f", condition_result$test$fstat[1], condition_result$test$pvalue[1]),
    sprintf("F=%.2f, p=%.3f", interaction_result$test$fstat[1], interaction_result$test$pvalue[1]),
    sprintf("%.4f %s units", metrics_CT$cumulative_travel, value_col),
    sprintf("%.4f %s units", metrics_TS$cumulative_travel, value_col),
    sprintf("%.3f", metrics_TS$cumulative_travel / metrics_CT$cumulative_travel),
    sprintf("%.2f%%", metrics_CT$cv),
    sprintf("%.2f%%", metrics_TS$cv),
    sprintf("%.4f %s units", metrics_CT$range_change, value_col),
    sprintf("%.4f %s units", metrics_TS$range_change, value_col),
    sprintf("%.4f %s units", cDT, value_col),
    sprintf("%.4f %s units", avg_abs_diff, value_col),
    sprintf("%.4f %s units", max_divergence, value_col)
  ),
  Significance = c(
    ifelse(time_CT_result$test$pvalue[1] < 0.05, "Yes", "No"),
    ifelse(time_TS_result$test$pvalue[1] < 0.05, "Yes", "No"),
    ifelse(condition_result$test$pvalue[1] < 0.05, "Yes", "No"),
    ifelse(interaction_result$test$pvalue[1] < 0.05, "Yes", "No"),
    "", "", "", "", "", "", "", "", "", ""
  )
)

# write.csv(combined_summary, paste0(base_filename, "_combined_summary.csv"), row.names = FALSE)
# cat("\n✓ Combined summary saved to:", paste0(base_filename, "_combined_summary.csv"), "\n")

#---------------------------------------------------------------------------
# Table 2: Model Fit Summary
#---------------------------------------------------------------------------
model_fit_stats <- data.frame(
  Parameter = c(
    "R-squared",
    "Adjusted R-squared",
    "Residual Standard Error",
    "F-statistic (overall model)",
    "Model p-value",
    "Number of observations"
  ),
  Value = c(
    round(summary(fit_spline)$r.squared, 4),
    round(summary(fit_spline)$adj.r.squared, 4),
    round(summary(fit_spline)$sigma, 4),
    round(summary(fit_spline)$fstatistic[1], 2),
    format(pf(summary(fit_spline)$fstatistic[1], 
              summary(fit_spline)$fstatistic[2], 
              summary(fit_spline)$fstatistic[3], 
              lower.tail = FALSE), scientific = TRUE, digits = 3),
    nrow(gi_summary)
  )
)

# write.csv(model_fit_stats, paste0(base_filename, "_model_fit.csv"), row.names = FALSE)
# cat("✓ Model fit statistics saved to:", paste0(base_filename, "_model_fit.csv"), "\n")

# cat("\nResults saved to directory:", output_dir, "\n")

#===============================================================================
# 8. VISUALIZATION
#===============================================================================

# Prepare prediction data
pred_df <- expand.grid(
  tp = seq(min(gi_summary$tp), max(gi_summary$tp), length.out = 200),
  condition = levels(gi_summary$condition)
)
pred_df[[value_col]] <- predict(fit_spline, newdata = pred_df)

plot_data <- data.frame(
  tp = tp_grid,
  CT = pred_CT,
  TS = pred_TS,
  Difference = diff_trajectory
)

# Main trajectory plot
p_main <- ggplot(gi_summary, aes(x = tp, y = .data[[value_col]], color = condition)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  geom_point(aes(shape = experiment), alpha = 0.5, size = 2) +
  geom_line(data = pred_df, aes(x = tp, y = .data[[value_col]], color = condition), size = 1.2) +
  scale_color_manual(values = color_mapping_condition, name = "Condition") +
  scale_shape_manual(values = 1:nlevels(gi_summary$experiment), name = "Experiment") +
  labs(
    title = paste(display_name, "Index Over Time"),
    subtitle = sprintf("CT: cT = %.3f | TS: cT = %.3f | Ratio TS/CT = %.2f",
                       metrics_CT$cumulative_travel, metrics_TS$cumulative_travel,
                       metrics_TS$cumulative_travel / metrics_CT$cumulative_travel),
    x = "Time [h]", 
    y = paste0(display_name, " Index [%]")
  ) +
  scale_linetype_manual(
    values = c("TS to 32 °C" = "dashed"),
    name = "Event"
  ) +
  theme_tidy() +
  theme(legend.position = "bottom")

plot(p_main)
# Difference plot
p_diff <- ggplot(plot_data, aes(x = tp, y = Difference)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  geom_line(color = "black", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Condition Difference (TS - CT)",
    subtitle = sprintf("cDT = %.3f | Avg diff = %.3f | Max divergence = %.3f",
                       cDT, avg_abs_diff, max_divergence),
    x = "Time [h]", 
    y = paste0("Difference in ", value_col)
  ) +
  scale_linetype_manual(
    values = c("TS to 32 °C" = "dashed"),
    name = "Event"
  ) +
  theme_tidy() +
  theme(legend.position = "")

# Display plots
p_combined <- grid.arrange(p_main, p_diff, ncol = 1)

# Save combined main and difference plot
# ggsave(file.path(output_fig_dir, paste0(index_type, "_fitted_trajectories.png")), 
       # p_combined, width = 8, height = 10, dpi = 300, bg = "white")
# ggsave(file.path(output_fig_dir, paste0(index_type, "_fitted_trajectories.pdf")), 
#        p_combined, width = 8, height = 10)
# 
# # Also save individual plots separately
# ggsave(file.path(output_fig_dir, paste0(index_type, "_main_trajectory.png")), 
#        p_main, width = 8, height = 6, dpi = 300,bg = "white")
# ggsave(file.path(output_fig_dir, paste0(index_type, "_main_trajectory.pdf")), 
#        p_main, width = 8, height = 6)
# 
# ggsave(file.path(output_fig_dir, paste0(index_type, "_difference_plot.png")), 
#        p_diff, width = 8, height = 4, dpi = 300,bg = "white")
# ggsave(file.path(output_fig_dir, paste0(index_type, "_difference_plot.pdf")), 
#        p_diff, width = 8, height = 4)

#===============================================================================
# 8. VISUALIZATION WITH AESTHETICS FROM METABOLITE PLOTS
#===============================================================================

# Define aesthetics from previous code
PLOT_FONT_SIZE <- 6
FONT_FAMILY <- "Liberation Sans"
LINE_WIDTH <- 0.6
AXIS_LINE_WIDTH <- 0.3

color_condition <- c(
  "CT" = "#D94801",
  "TS" = "#6A51A3"
)
event_color <- "#58A787"

# Temperature shift time (in hours)
temp_shift_h <- 146

# Prepare prediction data with same style as previous code
tp_grid <- seq(min(gi_summary$tp), max(gi_summary$tp), length.out = 200)
pred_CT <- predict(fit_spline, data.frame(tp = tp_grid, condition = "CT"))
pred_TS <- predict(fit_spline, data.frame(tp = tp_grid, condition = "TS"))

# Create prediction dataframe
pred_df <- data.frame(
  tp = rep(tp_grid, 2),
  value = c(pred_CT, pred_TS),
  condition = rep(c("CT", "TS"), each = length(tp_grid))
)

# Get significance flags from your test results
sig_flags <- data.frame(
  time_ct = time_CT_result$test$pvalue[1] < 0.05,
  time_ts = time_TS_result$test$pvalue[1] < 0.05,
  condition = condition_result$test$pvalue[1] < 0.05,
  interaction = interaction_result$test$pvalue[1] < 0.05
)

# Position for significance symbols
x_pos <- max(gi_summary$tp) * 0.85
y_base <- max(gi_summary[[value_col]]) * 1.05
y_step <- diff(range(gi_summary[[value_col]])) * 0.12

# Create significance symbols
sig_symbols <- data.frame(
  x = numeric(),
  y = numeric(),
  label = character(),
  color = character(),
  stringsAsFactors = FALSE
)

if(sig_flags$time_ct) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u25CF", 
                                  color = unname(color_condition["CT"])))
  y_base <- y_base - y_step
}
if(sig_flags$time_ts) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u25CF", 
                                  color = unname(color_condition["TS"])))
  y_base <- y_base - y_step
}
if(sig_flags$condition) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u0394", 
                                  color = "grey25"))
  y_base <- y_base - y_step
}
if(sig_flags$interaction) {
  sig_symbols <- rbind(sig_symbols, 
                       data.frame(x = x_pos, 
                                  y = y_base, 
                                  label = "\u00D7", 
                                  color = "grey25"))
}

# Main trajectory plot with metabolite aesthetics
p_main_aesthetic <- ggplot() +
  # # Raw data points
  # geom_point(data = gi_summary, 
  #            aes(x = tp, y = .data[[value_col]], color = condition),
  #            size = 1.5, alpha = 0.6) +
  # Spline lines
  geom_line(data = pred_df,
            aes(x = tp, y = value, color = condition),
            linewidth = LINE_WIDTH) +
  # Temperature shift line
  geom_vline(xintercept = temp_shift_h, 
             linetype = "dashed", 
             color = event_color, 
             linewidth = AXIS_LINE_WIDTH) +
  # Significance symbols
  geom_text(
    data = sig_symbols,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    color = sig_symbols$color,
    size = PLOT_FONT_SIZE / .pt,
    family = FONT_FAMILY,
    vjust = 0.5,
    hjust = 0.5
  ) +
  # Colors
  scale_color_manual(values = color_condition, name = "Condition") +
  # Labels
  labs(
    title = paste(display_name, "Index"),
    x = "Time (hours)",
    y = paste0(display_name, " Index [%]"),
    color = "Condition"
  ) +
  # Theme matching previous code
  theme_classic(base_size = PLOT_FONT_SIZE, base_family = FONT_FAMILY) +
  theme(
    plot.background = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = AXIS_LINE_WIDTH, colour = "black"),
    axis.ticks = element_line(linewidth = AXIS_LINE_WIDTH, colour = "black"),
    text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    plot.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY, 
                              face = "bold", hjust = 0.5),
    axis.text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    axis.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    legend.position = "bottom",
    legend.title = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY),
    legend.text = element_text(size = PLOT_FONT_SIZE, family = FONT_FAMILY)
  )

# Display the plot
print(p_main_aesthetic)


# Save main plot
ggsave(
  file.path(output_fig_dir, paste0(index_type, "_main_aesthetic.png")),
  plot = p_main_aesthetic,
  width = 80,
  height = 60,
  units = "mm",
  dpi = 300
)



#===============================================================================
# 9. MODEL DIAGNOSTICS
#===============================================================================

# Residual diagnostics
gi_summary$resid <- resid(fit_spline)
gi_summary$fitted <- fitted(fit_spline)

p_resid <- ggplot(gi_summary, aes(fitted, resid, color = condition)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = color_mapping_condition, name = "Condition") +
  labs(
    title = "Residual Plot",
    x = "Fitted values", 
    y = "Residuals"
  ) +
  theme_minimal()

# Q-Q plot
p_qq <- ggplot(gi_summary, aes(sample = resid)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line() +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles", 
    y = "Sample Quantiles"
  ) +
  theme_minimal()

# Display diagnostic plots
p_diag_combined <- grid.arrange(p_resid, p_qq, ncol = 2)

# Save combined diagnostic plots
# ggsave(file.path(output_fig_dir, paste0(index_type, "_diagnostics.png")), 
#        p_diag_combined, width = 10, height = 5, dpi = 300, bg = "white")
# ggsave(file.path(output_fig_dir, paste0(index_type, "_diagnostics.pdf")), 
#        p_diag_combined, width = 10, height = 5)



cat("\n========== ANALYSIS COMPLETE ==========\n")
