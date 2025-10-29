library(tidyverse)

octet_data <- read_csv("data/20240604_CharRun_Sum_ViCell_Titer_R_export.csv") %>%
  select(TP, Experiment, Hours, Titer, Group) %>%
  filter(Titer != "#N/A") %>%
  mutate(across(Titer, as.numeric)) %>%
  mutate(TP = str_replace(TP, "TP", "")) %>%
  mutate(across(TP, as.numeric))

# Define the custom x-axis labels
x_labels <- c(1, 8, 24, 32, 48, 56, 72, 74, 76, 80, 96, 100, 104, 120, 122, 124, 
              128, 144, 146, 148, 152, 168, 170, 172, 176, 192, 196, 200, 216, 
              218, 220, 224, 240, 244, 248, 264, 266, 268, 272, 288, 292, 296, 
              312, 314, 316, 320, 336, 340, 344, 360, 362, 364, 368, 384, 392, 
              408)

ggplot(data = octet_data, mapping = aes(x = Hours, y = Titer, group = Experiment)) +
  geom_line(aes(color = Experiment)) +
  geom_point(aes(color = Experiment)) +
  xlab("Time [hours]") +
  ylab("Titer [µg/ml]") +
  scale_x_continuous(
    name = "Time [hours]",
    breaks = x_labels,
    labels = x_labels,
    minor_breaks = seq(0, max(octet_data$Hours, na.rm = TRUE), by = 8)
  )  +
  scale_y_continuous(
    name = "Titer [µg/ml]",
    breaks = seq(0, max(octet_data$Titer, na.rm = TRUE), by = 100),
    minor_breaks = seq(0, max(octet_data$Titer, na.rm = TRUE), by = 100)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) 

unique(round(octet_data$Hours, digits = 0))

ggsave(filename = "figures/titer.png",    
       height = 160,
       width = 500,
       units = "mm",
       dpi = 600)

selected_titer <- octet_data %>%
  filter(TP %in% c(14, 28, 39, 46)) %>%

