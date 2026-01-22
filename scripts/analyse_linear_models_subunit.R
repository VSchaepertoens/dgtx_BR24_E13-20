library(tidyverse)
library(splines)
library(here)
source("scripts/plotting_utils.R")

input_file_path <- here::here("analysis", "charrun_E13-E20_subunit_V02_20260121_VS.RData")

load(file = input_file_path)

input_file_path <- here::here("analysis", "charrun_E13-E20_subunit_V01_20260121_VS.RData")

load(file = input_file_path)


# make df from the matrix -------------------------------------------------

clr_data_summarized <- clr_data.matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "subunit") %>%
  pivot_longer(cols = -subunit,
               names_to = "experiment_tp",
               values_to = "clr_fractional_abundance") %>%
  separate(col = experiment_tp,
           into = c('experiment', 'timepoint'),
           sep = "_",
           remove = FALSE
  ) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
    TRUE ~ 'other'
  ),
  timepoint = as.numeric(timepoint)) 

#  visualise data-----------------------------------------------------------

## 1. check linearity
## visualise whether relationship between response and explanatory variables is linear
ggplot(clr_data_summarized, aes(timepoint, clr_fractional_abundance, color = condition)) +
  geom_point(aes(shape = experiment)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed") +
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition),
    name = "Condition"
  ) +
  facet_wrap(~subunit) +
  scale_shape_manual(values = 1:nlevels(clr_data_summarized$experiment))

  
  
  filt_data <- clr_data_summarized %>% filter(subunit == "Intact") 
  
  filt_data$condition <- factor(filt_data$condition)
  levels(filt_data$condition)
  filt_data$condition <- relevel(filt_data$condition, ref = "CT")
  levels(filt_data$condition)
  
# formal test for linearity
fit_lin <- lm(clr_fractional_abundance ~ timepoint * condition, data = filt_data) # fit linear model
summary(fit_lin)
fit_quad <- lm(clr_fractional_abundance ~ (timepoint + I(timepoint^2)) * condition, data = filt_data) # fit quadratic model
summary(fit_quad )

anova(fit_lin, fit_quad) # test whether the quadratic term improves the fit

# formal test using splines
fit_spline <- lm(clr_fractional_abundance ~ ns(timepoint, df = 3) * condition, data = filt_data)
summary(fit_spline)

anova(fit_lin, fit_spline) # test whether the splines improve the fit

## 2. check if outliers present
ggplot(filt_data, aes(y = clr_fractional_abundance, x = condition)) +
  geom_boxplot()

# perform linear regression using splines ---------------------------------------------

# test for interaction of time and the condition
fit_spline <- lm(clr_fractional_abundance ~ ns(timepoint, df = 3) * condition, data = filt_data)
summary(fit_spline)

# test for no interaction of time and condition
fit_no_int <- lm(clr_fractional_abundance ~ ns(timepoint, df = 3) + condition, data = filt_data)
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
  timepoint = seq(min(filt_data$timepoint),
              max(filt_data$timepoint),
              length.out = 200),
  condition = levels(filt_data$condition)
)

# get fitted values
pred_df$spline_fit <- predict(fit_spline, newdata = pred_df)

filt_data$experiment <- factor(filt_data$experiment)

ggplot(filt_data, aes(timepoint, clr_fractional_abundance, color = condition)) +
  geom_vline(aes(xintercept = 146, linetype = "Temp. shift"),
             color = "#58A787", linewidth = 1) +
  geom_point(data = filt_data,
             aes(shape = experiment),
             alpha = 0.5) +
  geom_line(data = pred_df, aes(timepoint, spline_fit), linewidth = 1) +
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition),
    name = "Fitted splines per condition"
  ) +
  scale_shape_manual(values = 1:nlevels(filt_data$experiment),
                     name = "Experiment") +
  labs(
    x = "Time [h]",
    y = "CLR fractional abundance - Intact",
    linetype = "Temp. shift"
  )


# create residual plots, check that residuals are homoscedastic------------------------------------

par(mfrow = c(2, 2))
plot(fit_spline)
par(mfrow = c(1, 1))

filt_data$resid <- resid(fit_spline)
filt_data$fitted <- fitted(fit_spline)

ggplot(filt_data, aes(fitted, resid, color = condition)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals")













