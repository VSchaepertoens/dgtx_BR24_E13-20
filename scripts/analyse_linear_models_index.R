library(tidyverse)
library(splines)
library(here)

input_file_path <- here::here("analysis", "charrun_E13-E20_glycationindex_V01_20260114_VS.RData")

load(file = input_file_path)
# load("analysis/charrun_E13-E20_galactosylationindex_V01_20260114_VS.RData")
# load("analysis/charrun_E13-E20_glycationindex_V01_20260114_VS.RData")

gi_summary

# set color scheme --------------------------------------------------------

color_mapping_condition <- c(
  "CT" = "#E6641E",
  "TS" = "#4B288C"
)

#  visualise data-----------------------------------------------------------

## 1. check linearity
## visualise whether relationship between response and explanatory variables is linear
ggplot(gi_summary, aes(tp, GI, color = condition)) +
  geom_point(shape = gi_summary$experiment) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition),
    name = "Condition"
  ) +
  scale_shape_manual(values = 1:nlevels(gi_summary$experiment))
  

# formal test for linearity
fit_lin <- lm(GI ~ tp * condition, data = gi_summary) # fit linear model
summary(fit_lin)
fit_quad <- lm(GI ~ (tp + I(tp^2)) * condition, data = gi_summary) # fit quadratic model
summary(fit_quad )

anova(fit_lin, fit_quad) # test whether the quadratic term improves the fit

# formal test using splines
fit_spline <- lm(GI ~ ns(tp, df = 3) * condition, data = gi_summary)
summary(fit_spline)

anova(fit_lin, fit_spline) # test whether the splines improve the fit

## 2. check if outliers present
ggplot(gi_summary, aes(y = GI, x = condition)) +
  geom_boxplot()

# perform linear regression using splines ---------------------------------------------

# check reference level to be constant
gi_summary$condition <- factor(gi_summary$condition)
levels(gi_summary$condition)
gi_summary$condition <- relevel(gi_summary$condition, ref = "CT")
levels(gi_summary$condition)

# test for interaction of time and the condition
fit_spline <- lm(GI ~ ns(tp, df = 3) * condition, data = gi_summary)
summary(fit_spline)

# test for no interaction of time and condition
fit_no_int <- lm(GI ~ ns(tp, df = 3) + condition, data = gi_summary)
summary(fit_no_int)

anova_result <- anova(fit_no_int,fit_spline)
anova_result
# extract results
# extract results
interaction_stats <- anova_result[2, ]

# extract values
F_value <- interaction_stats$F
df_num  <- interaction_stats$Df        # numerator df
df_den  <- interaction_stats$Res.Df    # denominator df
p_value <- interaction_stats$`Pr(>F)`

# format p-value
p_text <- ifelse(p_value < 0.001, "< 0.001", sprintf("= %.3g", p_value))

# conditional text
significance_text <- if (p_value < 0.05) "was significant" else "was not significant"

interpretation_text <- if (p_value < 0.05) {
  "indicating that the dynamics over time in the Temp. shifted condition differ significantly from those in the Constant condition."
} else {
  "indicating that the dynamics over time in the Temp. shifted condition do not differ significantly from those in the Constant condition."
}

# construct and display sentence
cat(sprintf(
  "The time × condition interaction %s (F(%d, %d) = %.2f, p %s), %s\n",
  significance_text, df_num, df_den, F_value, p_text, interpretation_text
))


# visualise the results with the fitted splines --------------------------
# create new data for prediction
pred_df <- expand.grid(
  tp = seq(min(gi_summary$tp),
           max(gi_summary$tp),
           length.out = 200),
  condition = levels(gi_summary$condition)
)

# get fitted values
pred_df$GI_fit <- predict(fit_spline, newdata = pred_df)

gi_summary$experiment <- factor(gi_summary$experiment)

ggplot(gi_summary, aes(tp, GI, color = condition)) +
  geom_vline(aes(xintercept = 146, linetype = "Temp. shift"),
             color = "#58A787", linewidth = 1) +
  geom_point(data = gi_summary,
             aes(shape = experiment),
             alpha = 0.5) +
  geom_line(data = pred_df, aes(tp, GI_fit), linewidth = 1) +
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition),
    name = "Fitted splines per condition"
  ) +
  scale_shape_manual(values = 1:nlevels(gi_summary$experiment),
                     name = "Experiment") +
  labs(
    x = "Time [h]",
    y = "Glycation index [%]",
    linetype = "Temp. shift"
  )


# create residual plots, check that residuals are homoscedastic------------------------------------

par(mfrow = c(2, 2))
plot(fit_spline)
par(mfrow = c(1, 1))

gi_summary$resid <- resid(fit_spline)
gi_summary$fitted <- fitted(fit_spline)

ggplot(gi_summary, aes(fitted, resid, color = condition)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals")













