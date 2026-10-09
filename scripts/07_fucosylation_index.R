library(ggnewscale)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(purrr)
library(broom)
library(ggnewscale)
source("scripts/plotting_utils.R", local = FALSE)

# load abundance data -----------------------------------------------------

load("analysis/charrun_E13-E20_CQA_V06_20261006_VS.RData")

# calculate GI ---------------------------------------------------------

# Function to count fucose from glycan names
count_fucose <- function(glycan) {
  case_when(
    grepl("F", glycan) ~ 1,
    TRUE ~ 0
  )
}

# Split into two glycans and count galactose residues
corr_abundance_data <- corr_abundance_data %>%
  mutate(
    f1 = sub("/.*", "", glycoform1),
    f2 = sub(".*/", "", glycoform1),
    f1_fuc = sapply(f1, count_fucose),
    f2_fuc = sapply(f2, count_fucose),
    denominator_f1 = case_when(
      grepl("none", f1) ~ 0,
      TRUE ~ 1
    ),
    denominator_f2 = 1,
    total_fuc = (f1_fuc + f2_fuc) * corr_abundance,
    total_sites = (denominator_f1 + denominator_f2) * corr_abundance
  )

# Calculate GI per condition_br_tp
gi_summary <- corr_abundance_data %>%
  filter(!experiment_tp %in% c("E17_288", "E17_312")) %>%
  group_by(experiment_tp) %>%
  summarise(
    total_fuc = sum(total_fuc, na.rm = TRUE),
    total_sites = sum(total_sites, na.rm = TRUE)
  ) %>%
  mutate(
    GI = total_fuc / total_sites * 100
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


save(gi_summary, gi_stats, file = "analysis/charrun_E13-E20_fucosylationindex_V03_20261006_VS.RData")
load("analysis/charrun_E13-E20_fucosylationindex_V01_20260114_VS.RData")

# plot as dotplots over time --------------------------------------------------------

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
      # color = experiment
    ),
    linewidth = 0.6,
    alpha = 0.25
  ) +
  geom_point(
    data = gi_summary,
    aes(
      x = tp,
      y = GI,
      # color = experiment,
      shape = experiment
    ),
    size = 2,
    alpha = 0.25
  ) +
  # scale_color_manual(
  #   values = color_mapping_experiment,
  #   breaks = names(color_mapping_experiment),
  #   name = "Experiment"
  # ) +
  labs(
    color = "",
  ) +
  # new_scale_color() +
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
    y = "Fucosylation index [%]",
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


ggsave(filename = "figures/fucosylation_index_lineplot.png",
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
          file = "analysis/charrun_E13-E20_fucosylationindex_V01_20260114_VS.csv")

