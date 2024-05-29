library(tidyverse)
library("ggradar")

aa_rates <- read_csv("data/AA_rates.csv")

aa_rates <- aa_rates %>%
  select(-matches("SE"))

plot_radarplot <- function(data,
                           tshift = TRUE) {
  data_to_plot <- data %>%
    select(-1) %>%
    rename(group = Window) %>%
    column_to_rownames(var = "group") %>%
    select(all_of(sort(names(data[,3:20])))) %>%
    rownames_to_column(var = "group") %>%
    mutate(group = factor(group, 
                          levels = c("76 - 120 h",
                                      "124 - 168 h",
                                      "172 - 216 h",
                                      "220 - 264 h",
                                      "268 - 312 h")))
  if (tshift) {
    plot_title <- c("Temperature shift")
    color_scheme = c("#feb24c", "#fec44f", "#fe9929", "#ec7014", "#cc4c02")
    } else {
      plot_title <- c("No temperature shift")
      color_scheme = c("#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443")
      }
    
  ggradar(
    data_to_plot, 
    values.radar = c("-0.15", "0", "0.1"),
    grid.min = -0.15, grid.mid = 0, grid.max = 0.1,
    group.line.width = 1, 
    group.point.size = 3,
    group.colours = color_scheme,
    background.circle.colour = "white",
    gridline.mid.colour = "grey"
    ) +
    ggtitle(plot_title) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10) 
    )
}

plot_radarplot(aa_rates[6:10,])
ggsave("figures/aa_radarplot_tshifted.png",
       units = c("cm"),
       height = 15,
       width = 20,
       dpi = 600)

plot_radarplot(aa_rates[1:5,], tshift = FALSE)
ggsave("figures/aa_radarplot_nottshifted.png",
       units = c("cm"),
       height = 15,
       width = 20,
       dpi = 600)
