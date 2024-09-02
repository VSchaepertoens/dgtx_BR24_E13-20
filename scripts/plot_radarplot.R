library(tidyverse)
library("ggradar")
library(RColorBrewer)

brewer.pal(n = 9,"Oranges")
brewer.pal(n = 9,"Purples")


aa_rates <- read_csv("data/AA_rates.csv")

aa_rates <- aa_rates %>%
  select(-matches("SE"))

plot_radarplot <- function(data,
                           tshift = TRUE,
                           ncol = 20,
                           aa = TRUE) {
  data_to_plot <- data %>%
    select(-1) %>%
    rename(group = Window) %>%
    column_to_rownames(var = "group") %>%
    select(all_of(sort(names(data[,3:ncol])))) %>%
    rownames_to_column(var = "group") %>%
    mutate(group = factor(group, 
                          levels = c("76 - 120 h",
                                      "124 - 168 h",
                                      "172 - 216 h",
                                      "220 - 264 h",
                                      "268 - 312 h")))
  if (tshift) {
    plot_title <- c("Temperature shift")
    color_scheme = c("#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603")
    } else {
      plot_title <- c("No temperature shift")
      color_scheme = c("#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F")
    }
  
  if (aa) {
    radar_plot_labels <- c("-0.15", "0", "0.1")
    grid_min <- -0.15
    grid_mid <- 0
    grid_max <- 0.1
  } else {
    radar_plot_labels <- c("-2.6", "0", "2.3")
    grid_min <- -2.626331
    grid_mid <- 0
    grid_max <- 2.300745
  }
    
  ggradar(
    data_to_plot, 
    values.radar = radar_plot_labels,
    grid.min = grid_min, grid.mid = grid_mid, grid.max = grid_max,
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
    ) +
    guides(colour = guide_legend(ncol = 2)) 
    
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


# cedex data --------------------------------------------------------------

cedex_rates <- read_csv("data/CEDEX_rates.csv")

cedex_rates <- cedex_rates %>%
  select(-matches("SE"))

plot_radarplot(cedex_rates[6:10,],ncol = 5, aa = FALSE)
ggsave("figures/cedex_radarplot_nottshifted.png",
       units = c("cm"),
       height = 15,
       width = 20,
       dpi = 600)

plot_radarplot(cedex_rates[1:5,],tshift = FALSE, ncol = 5, aa = FALSE)
ggsave("figures/cedex_radarplot_tshifted.png",
       units = c("cm"),
       height = 15,
       width = 20,
       dpi = 600)





