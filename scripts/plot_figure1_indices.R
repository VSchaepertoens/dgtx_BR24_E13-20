library(tidyverse)
library(here)
source("scripts/plotting_utils.R", local = FALSE)

input_file_path <- here::here("analysis", "charrun_E13-E20_glycationindex_V01_20260114_VS.RData")
load(file = input_file_path)

glycation_index <- gi_stats %>%
  mutate(index = "Glycation")

input_file_path <- here::here("analysis", "charrun_E13-E20_galactosylationindex_V01_20260114_VS.RData")
load(file = input_file_path)

galactosylation_index <- gi_stats %>%
  mutate(index = "Galactosylation")


input_file_path <- here::here("analysis", "charrun_E13-E20_fucosylationindex_V01_20260114_VS.RData")
load(file = input_file_path)

fucosylation_index <- gi_stats %>%
  mutate(index = "Fucosylation")


indices_df <- bind_rows(glycation_index,galactosylation_index,fucosylation_index)



# plot all indices in one figure ------------------------------------------


ggplot(indices_df,aes(x = tp, group = condition)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  geom_ribbon(
    aes(
      ymin = mean_GI - se_GI,
      ymax = mean_GI + se_GI,
      fill = condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    aes(
      y = mean_GI,
      color = condition
    ),
    linewidth = 0.9,
    linetype = "solid"
  ) +
  facet_wrap(~index, ncol = 1, scales = "free_y") +
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
       y = "Index [%]") +
  theme_tidy() +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) + 
 theme(legend.position = "none")

ggsave(here("figures/figure1_indices.png"), width = 3, height = 3, dpi = 600, bg = "white")
ggsave(here("figures/figure1_indices.pdf"), width = 3, height = 3, dpi = 600, bg = "white")
