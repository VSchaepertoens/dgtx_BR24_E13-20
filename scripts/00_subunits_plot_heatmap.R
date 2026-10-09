library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(svglite)
library(compositions)
source("scripts/plotting_utils.R", local = FALSE)
library(here)


# load data 20251002 & preprocess -----------------------------------------

data <- read_csv('data/uv_subunits_TB_allTP.csv')

data_summarized <- data %>%
  group_by(experiment, timepoint) %>%
  summarise(lc_mean = mean(LC),
            lc2_mean = mean(LC2),
            intact_mean = mean(Intact)) %>%
  mutate(total = rowSums(across(c(lc_mean, lc2_mean, intact_mean)))) %>%
  mutate(LC = (lc_mean/total)*100,
         LC2 = (lc2_mean/total)*100,
         Intact = (intact_mean/total)*100) %>%
  select(experiment, timepoint ,LC ,LC2 ,Intact) %>%
  pivot_longer(cols = c("LC", "LC2", "Intact"),
               names_to = c("subunit"),
               values_to = "peak_area") %>%
  mutate(timepoint = as.numeric(timepoint),
         experiment = factor(experiment, levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20"))) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
    TRUE ~ 'other'
  )) %>%
  filter(timepoint != 72) %>%
  mutate(experiment_tp = paste(experiment, timepoint, sep = "_")) 


# Line plot
ggplot(data_summarized, aes(x = timepoint, 
                            y = peak_area, 
                            color = subunit, 
                            group = subunit)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~experiment, nrow = 2) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Time", y = "Relative peak area (%)",
       title = "Subunit composition over time",
       color = "Subunit")

# color_mapping_condition <- c(
#   # "E13" = "#FD8D3C",
#   # "E14" = "#9E9AC8",
#   # "E15" = "#F16913",
#   # "E16" = "#807DBA",
#   # "E17" = "#D94801",
#   # "E18" = "#6A51A3",
#   "CT" = "#E6641E",
#   "TS" = "#4B288C"
# )

# Line plot, facet per subunit
ggplot(data_summarized, aes(x = timepoint, 
                            y = peak_area, 
                            color = condition)) +
  geom_vline(aes(xintercept = 146, linetype = "Temp. shift to 32 °C"),
             color = "#58A787", 
             linewidth = 1.5) +
  geom_point(aes(shape = experiment),
             size = 1,
             alpha = 0.5) +
  geom_line(aes(group = experiment), alpha = 0.3) + #“Trend lines show locally weighted regression fits (LOESS) with no confidence interval (se = FALSE).”
  geom_smooth(size = 1.2, se = FALSE, alpha = 0.9) +
  
  scale_color_manual(values = color_mapping_condition) +
  scale_shape_manual(values = 1:nlevels(data_summarized$experiment)) +
  
  facet_wrap(~subunit, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
    plot.background = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    legend.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal",
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
  ) +
  labs(x = "Time [h]", y = "Relative peak area [%]",
       title = "Subunit composition over time",
       shape = "Experiment",
       color = "Condition",
       linetype = "Event")



ggsave("figures/20251001_TB_cNISTCHO_CharRUns_rem_allTP/subunit_quantification/subunit_line.pdf",
       width = 7,
       height = 6,
       dpi = 600,
       bg = "white")

ggsave("figures/20251001_TB_cNISTCHO_CharRUns_rem_allTP/subunit_quantification/subunit_line.svg",
       width = 7,
       height = 6,
       dpi = 600,
       bg = "white")


# average of biological replicates ----------------------------------------

averaged_data_summarized <- data_summarized %>%
  filter(!experiment_tp %in% c("E17_288", "E17_312")) %>%
  group_by(subunit, condition, timepoint) %>%
  summarise(
    mean_peak_area = mean(peak_area, na.rm = TRUE),
    se_peak_area = sd(peak_area, na.rm = TRUE)/sqrt(n()) 
  ) %>%
  ungroup() 

ggplot(averaged_data_summarized,aes(x = timepoint, group = condition)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  geom_ribbon(
    aes(
      ymin = mean_peak_area - se_peak_area,
      ymax = mean_peak_area + se_peak_area,
      fill = condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    aes(
      y = mean_peak_area,
      color = condition
    ),
    linewidth = 0.9,
    linetype = "solid"
  ) +
  facet_wrap(~subunit, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(96, 360),
    breaks = seq(96, 360, 48)
  ) +
  scale_fill_manual(values = color_mapping_condition, name = "Condition") +
  scale_color_manual(values = color_mapping_condition, guide = "none") +
  scale_linetype_manual(
    values = c("TS to 32 °C" = "dashed"),
    name = "Event"
  ) +
  labs(x = "Time [h]",
       y = "Relative peak area [%]") +
  theme_tidy() +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) + theme(legend.position = "none")

ggsave(here("figures/figure1_subunits_E17_288_312rem.png"), width = 4, height = 3, dpi = 600, bg = "white")
ggsave(here("figures/figure1_subunits_E17_288_312rem.pdf"), width = 4, height = 3, dpi = 600, bg = "white")


# build matrix from the subunit data --------------------------------------

data_summarized %>%
  ungroup() %>%
  dplyr::summarise(n = dplyr::n(), .by = c(subunit, experiment_tp)) %>%
  dplyr::filter(n > 1)

data.matrix <- data_summarized %>%
  ungroup() %>%
  mutate(experiment_tp = paste(experiment, timepoint, sep = "_")) %>%
  select(subunit, experiment_tp, peak_area) %>%
  pivot_wider(values_from = peak_area,
              names_from = experiment_tp) %>%
  column_to_rownames('subunit') %>%
  as.matrix() 

save(data_summarized, data.matrix, file = "analysis/charrun_E13-E20_subunit_V01_20260121_VS.RData")


# transform to clr space --------------------------------------------------

meta <- tibble(sample_name = colnames(data.matrix)) %>%
  separate(col = sample_name,
           into = c('experiment', 'timepoint'),
           sep = "_",
           remove = FALSE
  ) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
    TRUE ~ 'other'  # This handles any other experiments, if applicable
  ))

# clr transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix

save(clr_data.matrix, meta, file = "analysis/charrun_E13-E20_subunit_V02_20260121_VS.RData")

# build matrix from the subunit data, filter E17_288 and E17_312 and transform to clr space-----------

data.matrix <- data_summarized %>%
  filter(!experiment_tp %in% c("E17_288", "E17_312")) %>%
  ungroup() %>%
  mutate(experiment_tp = paste(experiment, timepoint, sep = "_")) %>%
  select(subunit, experiment_tp, peak_area) %>%
  pivot_wider(values_from = peak_area,
              names_from = experiment_tp) %>%
  column_to_rownames('subunit') %>%
  as.matrix() 

meta <- tibble(sample_name = colnames(data.matrix)) %>%
  separate(col = sample_name,
           into = c('experiment', 'timepoint'),
           sep = "_",
           remove = FALSE
  ) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
    TRUE ~ 'other'  # This handles any other experiments, if applicable
  ))

# clr transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix

save(data.matrix,clr_data.matrix, meta, file = "analysis/charrun_E13-E20_subunit_V04_20261008_VS.RData")
