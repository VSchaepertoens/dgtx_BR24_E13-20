library(tidyverse)

octet_data <- read_csv("data/20240604_CharRun_Sum_ViCell_Titer_R_export.csv") %>%
  select(TP, Experiment, Hours, Titer, Group) %>%
  filter(Titer != "#N/A") %>%
  mutate(across(Titer, as.numeric)) %>%
  mutate(TP = str_replace(TP, "TP", "")) %>%
  mutate(across(TP, as.numeric))

ggplot(data = octet_data, mapping = aes(x = Hours, y = Titer, group = Experiment)) +
  geom_line(aes(color = Experiment)) +
  geom_point(aes(color = Experiment)) +
  xlab("Time [hours]") +
  ylab("Titer [µg/ml]") +
  theme_bw() +
  theme(#text = element_text(size = 20),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

ggsave(filename = "figures/titer_axes14.png",    
       height = 160,
       width = 200,
       units = "mm",
       dpi = 600)