library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# load_data ---------------------------------------------------------------

data <- read_csv('data/Subunit_quantificaction_ E13_20.csv')


# heatmap -----------------------------------------------------------------
## plot heatmap of raw data ------------------------------------------------

data.matrix <- as.matrix(data)

Heatmap(data.matrix)
Heatmap(data.matrix, col = plasma(100))
Heatmap(data.matrix, col = rev(rainbow(10)))

## calculate z-score & plot heatmap -------------------------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])

min(scaled.data.matrix)
max(scaled.data.matrix)
col_fun = colorRamp2(c(-1.3, 0, 1.6), c("green", "white", "red"))
Heatmap(scaled.data.matrix, col = col_fun)
Heatmap(scaled.data.matrix, col = plasma(100))
Heatmap(scaled.data.matrix, col = rev(rainbow(10)))


# build correct table -----------------------------------------------------

data <- data %>%
  select(Sample, `LC (means)`, `LC2 (means)`, `Intact mAb (means)`) %>%
  pivot_longer(cols = c("LC (means)", "LC2 (means)", "Intact mAb (means)"),
               names_to = c("subunit","means"),
               names_sep = " ",
               values_to = "peak_area") %>%
  select(!c("means"))

#convert 'subunit' to factor and specify level order ## important for plotting order
data$subunit <- factor(data$subunit, levels = c('Intact', 'LC2', 'LC'))

# plot stacked bar chart -------------------------------------------------
plot_subunits <-  function(data_to_plot,
                           experiment = "E" #default to plot all experimnets
                           ){
  
  data_to_plot <-  data_to_plot %>% filter(str_detect(Sample, experiment))
  
  ggplot(data_to_plot, aes(y = Sample, x = peak_area, fill = subunit)) + 
  geom_bar(stat = "identity", position = "fill", width = .7) +
  xlab("peak_area (%)") +
  ylab("CHO cell experiment") +
  scale_fill_brewer(palette = "Accent") +
  scale_y_discrete(limits = rev) +
  theme(text = element_text(size = 9, 
                            # face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        axis.text.y = element_text(margin = margin(r = 0)),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
  )
}

plot_subunits(data_to_plot = data)
ggsave("figures/subunit_quantification_minus17_all_experiments.png",        
       width = 8.89,
       height = 15,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "E13")
ggsave("figures/subunit_quantification_e13.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "E14")
ggsave("figures/subunit_quantification_e14.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "E15")
ggsave("figures/subunit_quantification_e15.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "E16")
ggsave("figures/subunit_quantification_e16.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

# plot_subunits(data_to_plot = data,
#               experiment = "E17")
# ggsave("figures/subunit_quantification_e17.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "E18")
ggsave("figures/subunit_quantification_e18.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "E19")
ggsave("figures/subunit_quantification_e19.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "E20")
ggsave("figures/subunit_quantification_e20.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "336")
ggsave("figures/subunit_quantification_336.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "288")
ggsave("figures/subunit_quantification_288.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "216")
ggsave("figures/subunit_quantification_216.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

plot_subunits(data_to_plot = data,
              experiment = "")
ggsave("figures/subunit_quantification_216.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)


# #plot stacked barchart with labels
# ggplot(data_averaged, aes(x = cell_variant, y = mean_peak_area, fill = subunit)) + 
#   geom_col(position = "fill") +
#   geom_text(aes(label = ifelse(mean_peak_area == 0, "", round(percent*100,0))),
#             position = position_fill(vjust = 0.5)) +
#   geom_hline(yintercept = 0, linewidth = .5) +
#   labs(y = "peak area (%)") +
#   labs(x = "CHO cell variant") +
#   scale_fill_brewer(palette = "Accent") +
#   coord_flip() +
#   theme_bw() +
#   theme(
#     legend.position = "top",
#     panel.grid.major.y = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.border = element_blank(),
#     axis.text.x = element_blank()
#     ) 
#  
# ggsave("figures/stacked_bar/stacked_bar_percent_all_SC.png", 
#        width = 8.89,
#        height = 6.5,
#        units = c("cm"),
#        dpi = 600)






















