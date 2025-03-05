library(tidyverse)

#dataset 1 -----------------------------------------------
aa_rates <- read_csv("data/AA_rates.csv") %>%
  mutate(Window = str_replace_all(Window, c("76 - 120 h" = "window_1",
                                            "124 - 168 h" = "window_2",
                                            "172 - 216 h" = "window_3",
                                            "220 - 264 h" = "window_4",
                                            "268 - 312 h" = "window_5"))
         )

aa_rates_reordered <- aa_rates %>%
  # mutate(across(Asn_SE, as.numeric)) %>%
  pivot_longer(cols = -c(Condition, Window), 
               names_to = "Variable", 
               values_to = "Value") %>%
  separate(Variable, into = c("amino_acid", "Type"), sep = "_", extra = "merge") %>%
  pivot_wider(names_from = Type, values_from = Value) %>%
  rename("qp" = `NA`, "se" = "SE", "condition" = "Condition", "window" = "Window")


aa_rates_reordered_tshifted <- aa_rates_reordered %>%
  filter(condition == "Temperature shifted")

aa_rates_reordered_nottshifted <- aa_rates_reordered %>%
  filter(condition == "Not temperature shifted")


write_csv(aa_rates_reordered, "data/aa_rates_reordered.csv")
write_csv(aa_rates_reordered_nottshifted, "data/aa_rates_reordered_nottshifted.csv")
write_csv(aa_rates_reordered_tshifted, "data/aa_rates_reordered_tshifted.csv")
  
#dataset 2 ----------------------------------------------------------------
rates <- read_csv("data/Rates_DCW_GrowthRate.csv")

rates_with_growth <- rates %>%
  distinct(Experiment, Window, Condition, Growth_rate.h, Growth_rate_SE.h) %>%
  rename(Rate_mM.gDCW.h = Growth_rate.h, SD_mM.gDCW.h = Growth_rate_SE.h) %>%
  mutate(AA_meta = "Growth_rate") %>%
  bind_rows(rates) %>%
  select(Experiment, Window, Condition, AA_meta, Rate_mM.gDCW.h, SD_mM.gDCW.h) %>%
  rename(Rate = Rate_mM.gDCW.h, SD = SD_mM.gDCW.h)

unique(rates_with_growth$Experiment)
unique(rates_with_growth$Condition)
unique(rates_with_growth$Window)
unique(rates_with_growth$AA_meta)

rates_with_growth_tshifted <- rates_with_growth %>%
  filter(Condition == "Temp. shifted")

rates_with_growth_nottshifted <- rates_with_growth %>%
  filter(Condition == "Constant")

write_csv(rates_with_growth, "data/aa_rates_reordered_data2.csv")
write_csv(rates_with_growth_nottshifted, "data/aa_rates_reordered_data2_nottshifted.csv")
write_csv(rates_with_growth_tshifted , "data/aa_rates_reordered_data2_tshifted.csv")


# plot rates within windows -----------------------------------------------
color_mapping_experiment <- c(
  "E13" = "#FD8D3C",
  "E14" = "#9E9AC8",
  "E15" = "#F16913",
  "E16" = "#807DBA",
  "E17" = "#D94801",
  "E18" = "#6A51A3",
  "E19" = "#A63603",
  "E20" = "#54278F"
)

rates_with_growth <- read_csv("data/rates_condition_specific/aa_rates_reordered_data2.csv")

# Ensure that 'Window' is a factor to correctly handle the x-axis and dodging
rates_with_growth$Window <- factor(rates_with_growth$Window)

# Create the plot
ggplot(data = rates_with_growth, aes(x = Window, y = Rate, color = Experiment)) +
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust points with dodge to separate experiments
  geom_errorbar(
    aes(
      ymin = Rate - SD,  # Lower bound of error bar
      ymax = Rate + SD   # Upper bound of error bar
    ),
    position = position_dodge(width = 0.5),  # Ensure error bars are dodged like points
    width = 0.2,  # Control the width of the error bars (adjust as necessary)
    linewidth = 0.25  # Control line thickness
  ) +
  scale_color_manual(values = color_mapping_experiment) +  # Customize colors
  facet_wrap(~ AA_meta, nrow = 6, scales = "free_y") + # Facet by AA_meta, independent y scales
  theme_minimal()

ggsave(filename = "figures/aa_metabolites/aa_rates_reordered_data2.png",
       width = 210,
       height = 150, 
       units = "mm",
       dpi = 600,
       bg = "white")
  
  # Calculate variance for each combination of Window and AA_meta (across all experiments)
  variance_per_window_aa <- rates_with_growth %>%
    group_by(Window, AA_meta) %>%
    summarise(variance = var(Rate, na.rm = TRUE), .groups = "drop")
  
  # View the variance data
  print(variance_per_window_aa)  

  # Create a separate plot for variance
  ggplot(variance_per_window_aa, aes(x = Window, y = variance, group = AA_meta)) +
    geom_point(size = 3) +  # Plot variance as points
    geom_line() +  # Connect points with a line
    facet_wrap(~ AA_meta, nrow = 6, scales = "free_y") +  # Facet by AA_meta with independent y scales
    # scale_color_manual(values = color_mapping_experiment) +  # Customize colors for AA_meta
    theme_minimal() +
    labs(y = "Variance", x = "Window")  # Label axes for clarity
  

# window_specific dataset -------------------------------------------------

rates <- read_csv("data/rates_05032025/rates_mM_gDCW_h_new.csv")
  
  # Create the plot
  ggplot(data = rates, aes(x = Window, y = rate_mM_gDCW_h, color = Condition)) +
    geom_point(position = position_dodge(width = 0.5)) +  # Adjust points with dodge to separate experiments
    geom_errorbar(
      aes(
        ymin = rate_mM_gDCW_h - SD_mM_gDCW_h,  # Lower bound of error bar
        ymax = rate_mM_gDCW_h + SD_mM_gDCW_h   # Upper bound of error bar
      ),
      position = position_dodge(width = 0.5),  # Ensure error bars are dodged like points
      width = 0.2,  # Control the width of the error bars (adjust as necessary)
      linewidth = 0.25  # Control line thickness
    ) +
    # scale_color_manual(values = color_mapping_experiment) +  # Customize colors
    facet_wrap(~ Metabolite, nrow = 6, scales = "free_y") + # Facet by AA_meta, independent y scales
    theme_minimal()
  
  ggsave(filename = "figures/aa_metabolites/rates_mM_gDCW_h_new.png",
         width = 210,
         height = 150, 
         units = "mm",
         dpi = 600,
         bg = "white")
