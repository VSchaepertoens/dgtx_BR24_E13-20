# calculation of galactosylation index for each technical replicate and then plotting the mean + sd
library(tidyverse)
library(ggpubr)
library(dplyr)
library(purrr)
library(broom)
source("scripts/plotting_utils.R", local = FALSE)

# load abundance data -----------------------------------------------------

load("analysis/charrun_E13-E20_CQA_V02_20251124_VS.RData")

# calculate GI ---------------------------------------------------------

# Function to count galactose from glycan names
count_galactose <- function(glycan) {
  case_when(
    grepl("G2", glycan) ~ 2,
    grepl("G1", glycan) ~ 1,
    TRUE ~ 0
  )
}

# Split into two glycans and count galactose residues
corr_abundance_data <- corr_abundance_data %>%
  mutate(
    g1 = sub("/.*", "", glycoform1),
    g2 = sub(".*/", "", glycoform1),
    g1_gal = sapply(g1, count_galactose),
    g2_gal = sapply(g2, count_galactose),
    denominator_g1 = case_when(
      grepl("none", g1) ~ 0,
      TRUE ~ 2
    ),
    denominator_g2 = 2,
    total_gal = (g1_gal + g2_gal) * corr_abundance,
    total_sites = (denominator_g1 + denominator_g2) * corr_abundance
  )

# Calculate GI per condition_br_tp
gi_summary <- corr_abundance_data %>%
  filter(!experiment_tp %in% c("E17_288", "E17_312")) %>%
  group_by(experiment_tp) %>%
  summarise(
    total_gal = sum(total_gal, na.rm = TRUE),
    total_sites = sum(total_sites, na.rm = TRUE)
  ) %>%
  mutate(
    GI = total_gal / total_sites * 100
  ) %>%
  separate(experiment_tp, sep = "_", into = c("experiment","tp"),remove = FALSE) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
    TRUE ~ 'other'
  ),
  tp = as.numeric(tp),
  experiment = factor(experiment, levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20")))


# Print the summary table
print(gi_summary)
# Compute mean and standard deviation

# Calculate summary stats per condition and timepoint
gi_stats <- gi_summary %>%
  group_by(condition, tp) %>%
  summarise(
    mean_GI = mean(GI),
    sd_GI = sd(GI),
    se_GI = sd(GI, na.rm = TRUE)/sqrt(n()), 
    .groups = "drop"
  ) 

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


save(gi_summary, gi_stats, file = "analysis/charrun_E13-E20_galactosylationindex_V02_20260602_VS.RData")
load("analysis/charrun_E13-E20_galactosylationindex_V01_20260114_VS.RData")


# plot gal index as barplots -------------------------------------------------

# ggplot(data = gi_stats, aes(x = tp, y = mean_GI)) +
#   geom_col(aes(fill = condition),
#            position = position_dodge(width = 0.9),
#            color = "black") +
#   geom_errorbar(
#     aes(
#       ymin = mean_GI - sd_GI,
#       ymax = mean_GI + sd_GI,
#       group = condition
#     ),
#     position = position_dodge(.9),
#     width = 5,
#     linewidth = .25
#   ) +
#   facet_wrap(~condition
#              ) +
#   scale_fill_manual(
#     values = color_mapping_condition,
#     breaks = levels(gi_stats$condition)
#   ) +
#   labs(
#     x = "Timepoint [h]",
#     y = "Galactosylation index [%]",
#     fill = "Strategy"
#   ) +
#   scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
#   theme_bw() +
#   theme(
#     text = element_text( 
#       size = 11,
#       family = "sans",
#       colour = "black"
#     ),
#     axis.line = element_line(),
#     axis.text = element_text(color = "black", size = 11),
#     axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
#     axis.title.y = element_text(hjust = 0.5, face = "bold"),
#     axis.title.x = element_text(hjust = 0.5, face = "bold"),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     panel.border = element_blank(),
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold"),
#     legend.text = element_text(),
#     legend.box = "horizontal"
#   ) +
#   NULL


# plot as dotplots over time --------------------------------------------------------
library(ggnewscale)

ggplot() +
  geom_vline(aes(xintercept = 146, linetype = "Temp. shift"),
             color = "#58A787", linewidth = 1) +
  # individual experiments
  geom_line(
    data = gi_summary,
    aes(
      x = tp,
      y = GI,
      group = experiment,
      color = experiment
    ),
    linewidth = 0.6,
    alpha = 0.25
  ) +
  geom_point(
    data = gi_summary,
    aes(
      x = tp,
      y = GI,
      color = experiment,
      shape = experiment
    ),
    size = 2,
    alpha = 0.25
  ) +
  scale_color_manual(
    values = color_mapping_experiment,
    breaks = names(color_mapping_experiment),
    name = "Experiment"
  ) +
  labs(
    color = "",
  ) +
  new_scale_color() +
  # condition means
  geom_line(
    data = gi_stats,
    aes(
      x = tp,
      y = mean_GI,
      group = condition,
      color = condition
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition),
    name = "Condition"
  ) +
  scale_shape_manual(values = 1:nlevels(gi_summary$experiment)) +
  labs(
    x = "Time [h]",
    y = "Galactosylation index [%]",
    color = "Condition",
    linetype = "Temp. shift",
    shape = "Experiment"
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.title.position = "top",
    legend.text = element_text(size = 10),
    legend.box = "horizontal"
  )
  
  
  ggsave(filename = "figures/galactosylation_index_lineplot.png",
         width = 170,
         height = 100,
         units = "mm",
         dpi = 600,
         bg = "white")

# make wider table --------------------------------------------------------
gi_stats_wider <- gi_stats %>% 
  mutate(
    across(starts_with("mean_GI"), ~ round(.x, 2)),
    across(starts_with("sd_GI"), ~ round(.x, 2))
  ) %>%
  select(condition,tp, mean_GI, sd_GI) %>%
  pivot_wider(values_from = c(mean_GI, sd_GI), 
                         names_from = tp,
                         names_glue = "{.value}_{tp}") %>%
  mutate(condition = factor(condition, levels = c("CT", "TS"))) %>%
  arrange(condition)

write_csv(gi_stats_wider,
          file = "analysis/charrun_E13-E20_galactosylationindex_V01_20260114_VS.csv")
  

