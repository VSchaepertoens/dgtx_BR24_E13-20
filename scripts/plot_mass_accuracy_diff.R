# script to calcualte and plot the mass accuracy
library(tidyverse)

data <- read_csv("data/mass_accuracy/Mass_accuracy_E13_E20_disulphides_rtable.csv")

peak_diff <- data %>%
  group_by(experiment, tp, glycoform, technical_replicate) %>%
  summarise(diff = mass_difference_within_peak[peak == "front"] - mass_difference_within_peak[peak == "tail"], 
            .groups = "drop") %>%
  mutate(experiment_tp = paste(experiment, tp, sep = "_"))

peak_diff_mean <- peak_diff %>%
  group_by(experiment, tp, glycoform) %>%
  summarise(mean_diff = mean(diff),
            sd_diff = sd(diff)) %>%
  mutate(experiment_tp = paste(experiment, tp, sep = "_"))
  

color_mapping_experiment_tp <- c(
  "E13_120" = "#FD8D3C",
  "E13_144" = "#FD8D3C",
  "E13_168" = "#FD8D3C",
  "E13_192" = "#FD8D3C",
  "E13_216" = "#FD8D3C",
  "E13_288" = "#FD8D3C",
  "E13_336" = "#FD8D3C",
  "E14_120" = "#9E9AC8",
  "E14_144" = "#9E9AC8",
  "E14_168" = "#9E9AC8",
  "E14_192" = "#9E9AC8",
  "E14_216" = "#9E9AC8",
  "E14_288" = "#9E9AC8",
  "E14_336" = "#9E9AC8",
  "E15_120" = "#F16913",
  "E15_144" = "#F16913",
  "E15_168" = "#F16913",
  "E15_192" = "#F16913",
  "E15_216" = "#F16913",
  "E15_288" = "#F16913",
  "E15_336" = "#F16913",
  # "E16" = "#807DBA",
  # "E17" = "#D94801",
  "E18_120" = "#6A51A3",
  "E18_144" = "#6A51A3",
  "E18_168" = "#6A51A3",
  "E18_192" = "#6A51A3",
  "E18_216" = "#6A51A3",
  "E18_288" = "#6A51A3",
  "E18_336" = "#6A51A3",
  "E19" = "#A63603",
  "E20" = "#54278F"
)

color_mapping_experiment <- c(
  "E13" = "#FD8D3C",
  "E14" = "#9E9AC8",
  "E15" = "#F16913",
  # "E16" = "#807DBA",
  # "E17" = "#D94801",
  "E18" = "#6A51A3",
  "E19" = "#A63603",
  "E20" = "#54278F"

)


# plot mean +/- sd ----------------------------------------------------------

ggplot(peak_diff_mean, aes(x = tp, y = mean_diff, color = experiment_tp, group = experiment)) +
  geom_hline(yintercept = 0, color = "black") +  # Highlight y = 0 line
  geom_vline(xintercept = 146, color = "#5EA38A", linetype = "dashed", linewidth = 1) +  # Highlight y = 0 line
  geom_point(aes(color = experiment),
             position = position_dodge(width = 10),  
             size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = mean_diff - sd_diff,
      ymax = mean_diff + sd_diff,
      group = experiment_tp,
      color = experiment
    ),
    position = position_dodge(width = 10), 
    width = 0.2, 
    linewidth = 0.5
  ) +
  scale_color_manual(values = color_mapping_experiment) +  # Customize colors
  scale_x_continuous(breaks = unique(peak_diff_mean$tp)) +  # Ensure ticks match tp values
  facet_wrap(~ glycoform, nrow = 3) +
  labs(
    title = "Mass of free disulfides",
    x = "Timepoint [h]",
    y = "Mass [Da]",
    color = "experiment"  # Rename legend
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),  # Axis labels size
    axis.text = element_text(size = 12),   # Axis tick labels size
    plot.title = element_text(size = 12),  # Plot title size
    legend.title = element_text(size = 12),  # Legend title size
    legend.text = element_text(size = 12),   # Legend text size
    legend.position = "top",
    legend.box = "horizontal",
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE)  # Force the legend to be in one row
  )


ggsave(filename = "figures/mass_deviations_peakfront_peaktail_points.png",
       width = 200,
       height = 200,
       units = "mm",
       dpi = 300,
       bg = "white")


ggplot(peak_diff_mean, aes(x = tp, y = mean_diff, color = experiment_tp, group = experiment)) +
  geom_hline(yintercept = 0, color = "black") +  # Highlight y = 0 line
  geom_vline(xintercept = 146, color = "#5EA38A", linetype = "dashed", linewidth = 1) +  # Vertical line
  geom_smooth(method = "lm", se = FALSE, aes(group = experiment, color = experiment), size = 0.5) +  # Trendline, increased size
  geom_point(aes(color = experiment),
             position = position_dodge(width = 1),  
             size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = mean_diff - sd_diff,
      ymax = mean_diff + sd_diff,
      group = experiment_tp,
      color = experiment
    ),
    position = position_dodge(width = 1), 
    width = 0.2, 
    linewidth = 0.5
  ) +
  scale_color_manual(values = color_mapping_experiment) +  # Customize colors
  scale_x_continuous(breaks = unique(peak_diff_mean$tp)) +  # Ensure ticks match tp values
  facet_wrap(~ glycoform, nrow = 3) +
  labs(
    title = "Mass of free disulfides",
    x = "Timepoint [h]",
    y = "Mass [Da]",
    color = "experiment"  # Rename legend
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),  # Axis labels size
    axis.text = element_text(size = 12),   # Axis tick labels size
    plot.title = element_text(size = 12),  # Plot title size
    legend.title = element_text(size = 12),  # Legend title size
    legend.text = element_text(size = 12),   # Legend text size
    legend.position = "top",
    legend.box = "horizontal",
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE)  # Force the legend to be in one row
  )

# Save the plot
ggsave(filename = "figures/mass_deviations_peakfront_peaktail_trendlines.png",
       width = 200,  # Width in mm
       height = 200,  # Height in mm
       units = "mm",  # Set units to mm
       dpi = 300,  # Set DPI for testing
       bg = "white")



