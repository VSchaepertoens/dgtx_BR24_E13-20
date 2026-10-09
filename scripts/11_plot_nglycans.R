source("scripts/plotting_utils.R", local = FALSE)
library(tidyverse)
library(here)


# loading the data --------------------------------------------------------

load("analysis/charrun_E13-E20_CQA_V05_20260728_VS.RData")



# changes to data, should be done earlier, in plot_corrected_abund --------

br_averaged_corr_abundance_data <- corr_abundance_data %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'CT',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'TS',
    TRUE ~ 'other'
  )) %>%
  group_by(glycoform1, condition, timepoint) %>%
  summarise(
    mean_frac_ab = mean(corr_abundance),
    se_frac_ab = sd(corr_abundance)/sqrt(n()) 
  ) %>%
  ungroup() %>%
  mutate(timepoint = as.numeric(as.character(timepoint)),
         glycoform1 = gsub("/", " · ", glycoform1)) %>%
  filter(glycoform1 %in% c("G0F · G0F", "G0 · G0F", "G0F · G1F", "G1F · G1F"))

# plot lineplots ----------------------------------------------------------

ggplot(br_averaged_corr_abundance_data,aes(x = timepoint, group = condition)) +
  geom_vline(
    aes(xintercept = 146, linetype = "TS to 32 °C"),
    color = "#58A787",
    linewidth = 1
  ) +
  geom_ribbon(
    aes(
      ymin = mean_frac_ab - se_frac_ab,
      ymax = mean_frac_ab + se_frac_ab,
      fill = condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    aes(
      y = mean_frac_ab,
      color = condition
    ),
    linewidth = 0.9,
    linetype = "solid"
  ) +
  facet_wrap(~glycoform1, ncol = 2, scales = "free_y") +
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
       y = "Fractional abundance [%]") +
  theme_tidy() +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) + theme(legend.position = "none")

ggsave(here("figures/figure1_nglycans.png"), width = 5, height = 3, dpi = 600, bg = "white")
ggsave(here("figures/figure1_nglycans.pdf"), width = 5, height = 3, dpi = 600, bg = "white")



