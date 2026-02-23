library(tidyverse)
library(splines)
library(here)
source(here::here("scripts", "plotting_utils.R"))

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
  scale_shape_manual(values = 1:length(clr_data_summarized$experiment))


  
  # subset data for a single feature
  filt_data <- clr_data_summarized %>% filter(subunit == "LC") 
  
  filt_data$condition <- factor(filt_data$condition)
  levels(filt_data$condition)
  filt_data$condition <- relevel(filt_data$condition, ref = "CT")
  levels(filt_data$condition)
  
# formal test for linearity
fit_lin <- lm(clr_fractional_abundance ~ timepoint * condition, data = filt_data) # fit linear model
summary(fit_lin)
fit_quad <- lm(filt_data$clr_fractional_abundance ~ (timepoint + I(timepoint^2)) * condition, data = filt_data) # fit quadratic model
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
# make sure condition is properly defined once
clr_data_summarized <- clr_data_summarized %>%
  mutate(
    condition = factor(condition),
    condition = relevel(condition, ref = "CT")
  )

# get all subunits automatically
subunits <- unique(clr_data_summarized$subunit)

# empty list to store results
results_list <- list()

for (s in subunits) {
  
  # subset data
  filt_data <- clr_data_summarized %>%
    filter(subunit == s)
  
  # fit models
  fit_spline <- lm(clr_fractional_abundance ~ ns(timepoint, df = 3) * condition,
                   data = filt_data)
  
  fit_no_int <- lm(clr_fractional_abundance ~ ns(timepoint, df = 3) + condition,
                   data = filt_data)
  
  # compare models
  anova_result <- anova(fit_no_int, fit_spline)
  
  # extract interaction row (second row)
  interaction_stats <- anova_result[2, ]
  
  # store results
  results_list[[s]] <- data.frame(
    subunit = s,
    df_num  = interaction_stats$Df,
    df_den  = interaction_stats$Res.Df,
    F_value = interaction_stats$F,
    p_value = interaction_stats$`Pr(>F)`
  )
}

# combine into single dataframe
results_df <- bind_rows(results_list)

# optional: adjust for multiple testing
results_df <- results_df %>%
  mutate(
    p_adj_BH   = p.adjust(p_value, method = "BH"),
    p_adj_holm = p.adjust(p_value, method = "holm")
  )

results_df

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
    y = "CLR fractional abundance - LC",
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













