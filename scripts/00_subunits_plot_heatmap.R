library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(svglite)
library(compositions)
source("scripts/plotting_utils.R", local = FALSE)
library(here)


# load_data ---------------------------------------------------------------

# data <- read_csv('data/Subunit_quantificaction_ E13_20.csv')


# load data 20251002 & preprocess -----------------------------------------

data <- read_csv('data/uv_subunits_TB.csv')

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
  


# heatmap -----------------------------------------------------------------
# ## plot heatmap of raw data 
# 
# data.matrix <- as.matrix(data)
# 
# Heatmap(data.matrix)
# Heatmap(data.matrix, col = plasma(100))
# Heatmap(data.matrix, col = rev(rainbow(10)))
# 
# ## calculate z-score & plot heatmap 
# 
# scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  
# 
# #check for sanity
# mean(data.matrix[1,])
# sd(data.matrix[1,])
# (data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
# (data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])
# 
# min(scaled.data.matrix)
# max(scaled.data.matrix)
# col_fun = colorRamp2(c(-1.3, 0, 1.6), c("green", "white", "red"))
# Heatmap(scaled.data.matrix, col = col_fun)
# Heatmap(scaled.data.matrix, col = plasma(100))
# Heatmap(scaled.data.matrix, col = rev(rainbow(10)))


# build correct table -----------------------------------------------------

# data <- data %>%
#   select(Sample, `LC (means)`, `LC2 (means)`, `Intact mAb (means)`) %>%
#   pivot_longer(cols = c("LC (means)", "LC2 (means)", "Intact mAb (means)"),
#                names_to = c("subunit","means"),
#                names_sep = " ",
#                values_to = "peak_area") %>%
#   select(!c("means"))
# 
# #convert 'subunit' to factor and specify level order ## important for plotting order
# data$subunit <- factor(data$subunit, levels = c('Intact', 'LC2', 'LC'))

# # plot stacked bar chart -------------------------------------------------
# plot_subunits <-  function(data_to_plot,
#                            experiment = "E", #default to plot all experimnets
#                            plot_condition = FALSE
#                            ){
#   
#   data_to_plot <-  data_to_plot %>% filter(str_detect(Sample, experiment))
#   
#   p <- ggplot(data_to_plot, aes(y = Sample, x = peak_area, fill = subunit)) + 
#       geom_bar(stat = "identity", position = "fill", width = .7) +
#       xlab("peak_area (%)") +
#       ylab("CHO cell experiment") +
#       scale_fill_brewer(palette = "Accent") +
#       scale_y_discrete(limits = rev) +
# 
#       theme(text = element_text(size = 16, 
#                                 face = "bold",
#                                 family = "sans"),
#             axis.text = element_text(colour = "black"),
#             panel.background = element_blank(),
#             axis.text.y = element_text(margin = margin(r = 0)),
#             axis.ticks.y = element_blank(),
#             legend.title = element_blank(),
#             legend.position = "bottom",
#             panel.border = element_blank(),
#             panel.grid.major.y = element_blank(),
#             panel.grid.minor = element_blank(),
#       )
#   
#   # Conditionally add the facet
#   if (plot_condition) {
#     p <- p + facet_wrap(~ condition)
#   }
#   # Return the plot
#   return(p)
#   
# }
# 
# plot_subunits(data_to_plot = data)
# ggsave("figures/subunit_quantification_minus17_all_experiments.png",        
#        width = 8.89,
#        height = 15,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "E13")
# ggsave("figures/subunit_quantification_e13.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "E14")
# ggsave("figures/subunit_quantification_e14.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "E15")
# ggsave("figures/subunit_quantification_e15.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "E16")
# ggsave("figures/subunit_quantification_e16.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# # plot_subunits(data_to_plot = data,
# #               experiment = "E17")
# # ggsave("figures/subunit_quantification_e17.png",        
# #        width = 8.89,
# #        height = 8.89,
# #        units = c("cm"),
# #        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "E18")
# ggsave("figures/subunit_quantification_e18.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "E19")
# ggsave("figures/subunit_quantification_e19.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "E20")
# ggsave("figures/subunit_quantification_e20.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "336")
# ggsave("figures/subunit_quantification_336.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "288")
# ggsave("figures/subunit_quantification_288.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# plot_subunits(data_to_plot = data,
#               experiment = "216")
# ggsave("figures/subunit_quantification_216.png",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)
# 
# 
# selected_data <- data %>%
#   filter(Sample %in% c("E13_120","E13_216", "E13_288", "E13_336",
#                        "E14_120","E14_216", "E14_288", "E14_336",
#                        "E19_120","E19_216", "E19_288", "E19_336",
#                        "E20_120","E20_216", "E20_288", "E20_336")) %>%
#   mutate(Sample = factor(Sample,
#                          levels = c("E13_120","E19_120","E14_120", "E20_120",  
#                                     "E13_216","E19_216","E14_216", "E20_216", 
#                                     "E13_288","E19_288","E14_288", "E20_288",
#                                     "E13_336","E19_336","E14_336", "E20_336"))) %>%
#   mutate(condition = case_when(
#     Sample %in% c("E13_120","E13_216", "E13_288", "E13_336","E19_120","E19_216", "E19_288", "E19_336") ~ 'constant',
#     Sample %in% c("E14_120","E14_216", "E14_288", "E14_336","E20_120","E20_216", "E20_288", "E20_336") ~ 'tshifted',
#     TRUE ~ 'other'  # This handles any other experiments, if applicable
#   ))
# 
# plot_subunits(data_to_plot = selected_data,
#               plot_condition = FALSE)
# ggsave("figures/subunit_quantification_4tp_4exp_tp_ordered.png",        
#        width = 15,
#        height = 15,
#        units = c("cm"),
#        dpi = 600)
# 
# 
# 
# 
# # #plot stacked barchart with labels
# # ggplot(data_averaged, aes(x = cell_variant, y = mean_peak_area, fill = subunit)) + 
# #   geom_col(position = "fill") +
# #   geom_text(aes(label = ifelse(mean_peak_area == 0, "", round(percent*100,0))),
# #             position = position_fill(vjust = 0.5)) +
# #   geom_hline(yintercept = 0, linewidth = .5) +
# #   labs(y = "peak area (%)") +
# #   labs(x = "CHO cell variant") +
# #   scale_fill_brewer(palette = "Accent") +
# #   coord_flip() +
# #   theme_bw() +
# #   theme(
# #     legend.position = "top",
# #     panel.grid.major.y = element_blank(),
# #     panel.grid.major.x = element_blank(),
# #     panel.grid.minor.x = element_blank(),
# #     axis.ticks.y = element_blank(),
# #     axis.ticks.x = element_blank(),
# #     panel.border = element_blank(),
# #     axis.text.x = element_blank()
# #     ) 
# #  
# # ggsave("figures/stacked_bar/stacked_bar_percent_all_SC.png", 
# #        width = 8.89,
# #        height = 6.5,
# #        units = c("cm"),
# #        dpi = 600)

# Line plots of relative abundances ---------------------------------------
# Assuming df has columns: Sample, subunit, peak_area
# all_conditions <- c("E17")

# df_clean <- data_summarized %>%
#   # separate sample into condition and time (assuming format E13_120)
#   # tidyr::separate(Sample, into = c("Experiment", "Time"), sep = "_") %>%
#   complete(Experiment = all_conditions,
#            Time,
#            subunit,
#            fill = list(rel_area = 0)) %>%
#   mutate(Time = as.numeric(Time),
#          Experiment = factor(Experiment, levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20"))) %>%
#   mutate(Condition = case_when(
#     Experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'Constant',
#     Experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'Temp. shifted',
#     TRUE ~ 'other'
#     )) 


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

# # Area plots (like stacked lines) -----------------------------------------
# 
# ggplot(df_clean, aes(x = Time, y = peak_area, fill = subunit)) +
#   geom_area(alpha = 0.8, color = "black", size = 0.2) +
#   facet_wrap(~Experiment, nrow = 2) +
#   theme_minimal(base_size = 11) +
#   theme(axis.text.x = element_text(angle = 45)) +
#   labs(x = "Time", y = "Relative peak area (%)",
#        title = "Stacked area plot of subunit composition",
#        fill = "Subunit")
# 
# # load titer data ---------------------------------------------------------
# 
# octet_data <- read_csv("data/20240604_CharRun_Sum_ViCell_Titer_R_export.csv") %>%
#   select(TP, Experiment, Hours, Titer, Group) %>%
#   filter(Titer != "#N/A") %>%
#   mutate(across(Titer, as.numeric)) %>%
#   mutate(TP = str_replace(TP, "TP", "")) %>%
#   mutate(across(TP, as.numeric))
# 
# selected_titer <- octet_data %>%
#   filter(TP %in% c(14, 28, 39, 46)) %>%
#   filter(Experiment %in% c("E13","E14", "E19", "E20")) %>%
#   mutate(hours_round = trunc(Hours)) %>%
#   mutate(experiment_tp = paste(Experiment, hours_round, sep = "_")) %>%
#   select(experiment_tp, Titer)
# 
# #need to use here a different timepoint, because the 336 timepoint is not present
# E13_selected_titer <- octet_data %>%
#   filter(TP %in% c(51)) %>%
#   filter(Experiment %in% c("E13")) %>%
#   mutate(hours_round = 336) %>%
#   # mutate(hours_round = trunc(Hours)) %>%
#   mutate(experiment_tp = paste(Experiment, hours_round, sep = "_")) %>%
#   select(experiment_tp, Titer)
# 
# selected_titer <- rbind(selected_titer,E13_selected_titer)
# 
# ggplot(selected_titer,aes(y = experiment_tp, x = Titer)) +
#   geom_point() +
#   geom_segment(aes(y = experiment_tp, yend = experiment_tp, x = 0, xend = Titer))
# 
# 
# 
# # combining the two plots -------------------------------------------------
# 
# library(ggplot2)
# library(dplyr)
# library(stringr)
# library(cowplot) # or use patchwork
# 
# # Your existing barplot
# p_bar <- ggplot(selected_data, aes(y = Sample, x = peak_area, fill = subunit)) + 
#   geom_bar(stat = "identity", position = "fill", width = 0.7) +
#   xlab("Proportion (%)") +
#   ylab("") +
#   scale_fill_brewer(palette = "Accent") +
#   scale_y_discrete(limits = rev(unique(selected_data$Sample))) + 
#   theme_minimal() +
#   theme(
#     text = element_text(size = 10, 
#                         face = "bold",
#                         family = "sans"),
#     axis.text = element_text(colour = "black"),
#     panel.background = element_blank(),
#     axis.text.y = element_text(margin = margin(r = 0)),
#     axis.ticks.y = element_blank(),
#     legend.title = element_blank(),
#     legend.position = "top",
#     panel.border = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
#  plot(p_bar)
# 
#  # Ensure the sample order matches
#  selected_titer$Sample <- factor(selected_titer$experiment_tp, levels = rev(unique(selected_data$Sample)))
#  
# 
#  p_lollipop <- ggplot(selected_titer, aes(x = Titer, y = Sample)) +
#    geom_segment(aes(xend = Titer, yend = Sample, x = 0), color = "black") +
#    # geom_linerange(aes(xmin = 0, xmax = Titer, y = Sample), color = "black", linewidth = 0.5)+
#    geom_point(size = 3, color = "black") +
#    scale_x_reverse(position = "top") + 
#    scale_y_discrete(limits = rev(unique(selected_data$Sample))) +
#    xlab("Titer (µg/mL)") +
#    ylab(NULL) +
#    theme_minimal() +
#    theme(
#      axis.text.y = element_blank(),  # Hide y-axis labels
#      # axis.ticks.y = element_blank(),
#      # panel.grid.major.y = element_blank(),
#      panel.grid.minor.x = element_blank(),
#      
#      text = element_text(size = 10, 
#                          face = "bold",
#                          family = "sans"),
#      axis.text = element_text(colour = "black"),
#      panel.background = element_blank(),
#      # axis.text.y = element_text(margin = margin(r = 0)),
#      axis.ticks.y = element_blank(),
#      legend.title = element_blank(),
#      legend.position = "bottom",
#      panel.border = element_blank(),
#      panel.grid.major.y = element_blank(),
#      panel.grid.minor = element_blank()
#    )
#  
# 
# 
#  plot(p_lollipop) 
# 
#  combined <- plot_grid(p_bar, 
#                        p_lollipop, 
#                        align = "h",
#                        nrow = 1, 
#                        rel_widths = c(1, 0.5))
#  
# 
#  print(combined)
# ggsave("figures/subunit_quantification/subunit_quantification_4tp_4exp_titer.png",
#        width = 6,
#        height = 4,
#        dpi = 300,
#        bg = "white")
# 
# ggsave("figures/subunit_quantification/subunit_quantification_4tp_4exp_titer.pdf",
#        width = 6,
#        height = 4,
#        dpi = 300,
#        bg = "white")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
