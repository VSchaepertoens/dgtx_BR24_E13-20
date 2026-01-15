library(ggnewscale)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(purrr)
library(broom)
library(ggnewscale)


# glycation index ---------------------------------------------------------
load("analysis/abundance_data_pngase.RData")

glycation_data <- abundance_data_averaged %>%
  filter(timepoint != "72")

# Check if frac_abundance sums to 1 for each replicate
abundance_sums <- glycation_data %>%
  group_by(experiment_tp) %>%
  summarise(total_abundance = sum(frac_abundance)) %>%
  ungroup()
#checked, all abundances sum to 100 

# Function to count glucose from glycation names
count_glucose <- function(glycation) {
  case_when(
    grepl("1xHex", glycation) ~ 1,
    grepl("2xHex", glycation) ~ 2,
    grepl("3xHex", glycation) ~ 3,
    TRUE ~ 0
  )
}

# Split into two glycans and count galactose residues
glycation_data <- glycation_data %>%
  mutate(
    glu = sapply(modcom_name, count_glucose),
    denominator_glu = 3,
    # denominator_glu = case_when(
    #   grepl("none", modcom_name) ~ 3,
    #   TRUE ~ 3),
    total_glu = glu * frac_abundance,
    # total_sites = denominator_glu * frac_abundance
  ) 
  # separate(experiment_tp, sep = "_", into = c("experiment","tp"), remove = FALSE) 

# sanity_check <- glycation_data %>%
#   group_by(experiment_tp) %>%
#   summarise(
#     max_possible_glu = sum(3 * frac_abundance),  # 3 is max glu per antibody
#     total_glu = sum(glu * frac_abundance)
#   ) %>%
#   mutate(percent = total_glu / max_possible_glu * 100) %>%
#   arrange(desc(percent))

#checked that all percent are below 100

# Calculate GI per condition_br_tp
gi_summary <- glycation_data %>%
  group_by(experiment_tp) %>%
  summarise(
    total_glu = sum(glu * frac_abundance, na.rm = TRUE),
    total_sites = sum(denominator_glu * frac_abundance, na.rm = TRUE)
  ) %>%
  mutate(
    glycation_index = (total_glu / total_sites) * 100
  )%>%
  separate(experiment_tp, sep = "_", into = c("experiment","tp"),remove = FALSE) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'Constant',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'Temp. shifted',
    TRUE ~ 'other'
  ),
  tp = as.numeric(tp),
  experiment = factor(experiment, levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20")))
  # mutate(time_group = if_else(tp == 120, "120", "240_264")) %>%
  # mutate(condition = case_when(
  #   condition == "A" ~ "STD",
  #   condition == "B" ~ "STD+",
  #   condition == "G" ~ "LoG",
  #   condition == "C" ~ "LoG+",
  #   condition == "D" ~ "HiF",
  #   condition == "E" ~ "HIP",
  #   condition == "F" ~ "HIP+")
  # ) 

# Print the summary table
print(gi_summary)
# Compute mean and standard deviation

# Calculate summary stats per condition and timepoint
gi_stats <- gi_summary %>%
  group_by(condition, tp) %>%
  summarise(
    mean_GI = mean(glycation_index),
    sd_GI = sd(glycation_index),
    .groups = "drop"
  ) 
  # mutate(time_group = if_else(tp == 120, "exponential", "stationary"),
  #        condition = factor(condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

save(gi_summary, gi_stats, file = "analysis/glycation_index.RData")
load("analysis/glycation_index.RData")

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
  mutate(condition = factor(condition, levels = c("Constant", "Temp. shifted"))) %>%
  arrange(condition)

write_csv(gi_stats_wider,
          file = "analysis/charrun_E13-E20_glycationnindex_V01_20260114_VS.csv")

# plot glycation index ----------------------------------------------------
color_mapping_condition <- c(
  "Constant" = "#A63603",
  "Temp. shifted" = "#54278F"
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

glu_ind <-
  ggplot(data = gi_stats) +
  # Points: color and linetype mapped to condition in one aes() call
  geom_point(
    aes(x = tp, y = mean_GI, color = condition),
    size = 2.5,
    position = position_dodge(width = 0.3)
  ) +
  
  # Lines: same combined mapping
  geom_line(
    aes(x = tp, y = mean_GI, color = condition),
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  
  # Error bars: also include linetype to ensure legend combines
  geom_errorbar(
    aes(
      x = tp,
      ymin = mean_GI - sd_GI,
      ymax = mean_GI + sd_GI,
      color = condition
      # linetype = condition
    ),
    width = 0.1,
    position = position_dodge(width = 0.3)
  ) +
  
  # Manual color and linetype mappings
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition)
  ) +
  
  # Unified legend title
  labs(
    x = "Time [h]",
    y = "Glycation index [%]",
    color = "Condition"
  ) +
  
  # scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  
  guides(color = guide_legend(nrow = 1)) +
  
  # Styling
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
    legend.text = element_text(size = 10),
    legend.box = "horizontal"
  )

plot(glu_ind)



# lineplots glycation index -----------------------------------------------
ggplot() +
  geom_vline(aes(xintercept = 146, linetype = "Temp. shift"),
             color = "#58A787", linewidth = 1) +
  # individual experiments
  geom_line(
    data = gi_summary,
    aes(
      x = tp,
      y = glycation_index,
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
      y = glycation_index,
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
    y = "Glycation index [%]",
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


ggsave(filename = "figures/glycation_index_lineplot.png",
       width = 170,
       height = 100,
       units = "mm",
       dpi = 600,
       bg = "white")


# # plot glycation index as barplot -----------------------------------------
# gly_ind_bar <- ggplot(data = gi_stats, aes(x = condition, y = mean_GI)) +
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
#     width = .5,
#     linewidth = .25
#   ) +
#   # geom_text(
#   #   aes(label = condition, fill = condition, y = 2),  # include fill here!
#   #   position = position_dodge(width = 0.9),
#   #   vjust = 0,
#   #   hjust = 0, 
#   #   angle = 90,
#   #   colour = "white",
#   #   size = 3
#   # ) +
#   facet_wrap(~time_group, 
#              labeller = labeller(time_group = c("exponential" = "Exponential",
#                                                 "stationary" = "Stationary"))
#   ) +
#   # # Lines: same combined mapping
#   # geom_line(
#   #   aes(x = time_group, y = mean_GI, color = condition, group = condition),
#   #   linewidth = 1,
#   #   position = position_dodge(width = 0.9)
#   # ) +
#   scale_fill_manual(
#     values = color_mapping_condition,
#     breaks = levels(gi_stats$condition)
#   ) +
#   # scale_color_manual(
#   #   values = color_mapping_condition,
#   #   breaks = names(color_mapping_condition)
#   # ) 
#   # Unified legend title
#   labs(
#     x = "",
#     y = "Glycation index (%)",
#     fill = "Strategy"
#   ) +
#   
#   scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
#   theme_bw() +
#   theme(
#     text = element_text( 
#       size = 11,
#       family = "sans",
#       colour = "black"
#     ),
#     axis.line = element_line(),
#     axis.text = element_text(color = "black", size = 11),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
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
#   
#   guides(fill = guide_legend(nrow = 1)) 
# 
# plot(gly_ind_bar)
# ggsave(filename = "figures/glycation_index_barplot_facet_time.png",
#        width = 150,
#        height = 100,
#        units = "mm",
#        dpi = 600,
#        bg = "white")
# 

