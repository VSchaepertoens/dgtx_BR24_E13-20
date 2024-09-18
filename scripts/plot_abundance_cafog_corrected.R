library(tidyverse)
#library("scales")    
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(fs)

# load cafog corrected data -----------------------------------------------

abundance_data <- NULL
# Specify the directory path
directory_path <- "analysis/cafog/"

# List all directories in the specified path
folders <- dir_ls(directory_path, type = "directory")

for (folder in folders) {
  file_path <- paste0(folder,"/results.csv")
  
  abundance_data <- rbind(abundance_data,
                          read_csv(file_path,
                                   n_max = 10) %>%
                            mutate(experiment_tp = str_extract(folder, "[^/]+$")) %>%
                            {.}
  )
}

corr_abundance_data <- abundance_data %>%
  separate(
    glycoform,
    into = c("glycoform1", "glycoform2",  "glycoform3"),
    sep = "\\s+or\\s+",
    remove = FALSE
  ) %>%
  select(glycoform1, corr_abundance, corr_abundance_error, experiment_tp) %>%
  mutate(glycoform1 = str_replace_all(glycoform1, "A2", ""))  %>%
  mutate(glycoform1 = str_replace_all(glycoform1, c("G0F/G2F" = "G1F/G1F", 
                                                    "G2F/none" = "none/G2F",
                                                    "G1F/none" = "none/G1F", 
                                                    "G0F/none" = "none/G0F", 
                                                    "G0/none" = "none/G0")
                                      )) %>%
  mutate(glycoform1 = factor(glycoform1, levels = c("none/G0F",
                                                    "none/G1F",
                                                    "none/G2F",
                                                    "G0/G0",
                                                    "G0/G0F",
                                                    "G0F/G0F",
                                                    "G0F/G1F",
                                                    "G1F/G1F",
                                                    "G1F/G2F",
                                                    "G2F/G2F")
                             )) %>%
  separate(experiment_tp,
           into = c("experiment",
                    "timepoint"
           ),
           sep = "_",
           remove = FALSE) %>%
  mutate(experiment = factor(experiment, levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20"))) %>%
  # filter(experiment != "E17")
  drop_na()
  {.}

save(corr_abundance_data, file = "analysis/corr_abundance_data.RData")

load("analysis/corr_abundance_data.RData")

corr_abundance_data <- corr_abundance_data %>%
  drop_na() %>%
  mutate(experiment = factor(experiment, 
                             levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20")
                             )) %>%
  filter(experiment != c("E17")) %>%
  filter(timepoint %in% c("120", "216", "288","336"))
# plot char runs data -----------------------------------------------------
# Define the colors
# color_mapping_experiment <- c(
#   "E13" = "#e41a1c",
#   "E14" = "#377eb8",
#   "E15" = "#4daf4a",
#   "E16" = "#984ea3",
#   "E17" = "#ff7f00",
#   "E18" = "#ffff33",
#   "E19" = "#a65628",
#   "E20" = "#f781bf"
# )

color_mapping_experiment <- c(
  "E13" = "#7f3b08",
  "E14" = "#2d004b",
  "E15" = "#b35806",
  "E16" = "#542788",
  "E17" = "#e08214",
  "E18" = "#8073ac",
  "E19" = "#fdb863",
  "E20" = "#b2abd2"
)


plot_bars <- function(data,
                      title = "Fractional abundance"){
  ggplot(data, aes(x = glycoform1, y = corr_abundance, fill = experiment)) +
    geom_col(
      position = position_dodge(width = 0.9)  
    ) + 
    geom_errorbar(
      aes(
        ymin = corr_abundance - corr_abundance_error,
        ymax = corr_abundance + corr_abundance_error,
        group = experiment
      ),
      position = position_dodge(.9),
      width = .5,
      linewidth = .25
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    facet_wrap(~ timepoint, nrow = 1) +
    scale_fill_manual(values = color_mapping_experiment) +
    scale_y_continuous(name = "fractional abundance (%)",
                       limits = c(0,65)) +
    xlab("") +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              face = "plain",
                              family = "sans"),
          axis.text.x = element_text(angle = 90, 
                                     vjust = .5, 
                                     hjust = 1),
          axis.text = element_text(colour = "black"),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          panel.border = element_blank()
          ) +
    ggtitle(title) +
    NULL
  }
  

#plot all experiments
plot_bars(corr_abundance_data,
          title = "Fractional abundance of glycans in all experiments")

ggsave(filename = "figures/corrected_frac_ab_barplot_minus17_all_experiments_4tp.png",
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)

#plot all experiments minus 17
# plot_bars(corr_abundance_data %>% filter(experiment != "E17"))
# 
# ggsave(filename = "figures/corrected_frac_ab_barplot_minus17_all_experiments.png",
#        height = 100,
#        width = 250,
#        units = "mm",
#        dpi = 600)

constant_data <- corr_abundance_data %>%
  filter(experiment %in% c("E13", "E15", "E19"))
plot_bars(constant_data,
          title = "Fractional abundance of glycans in constant temperature experiments")

ggsave(filename = "figures/corrected_frac_ab_barplot_minus17_constant_experiments_4tp.png",
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)

tshifted_data <- corr_abundance_data %>%
  filter(experiment %in% c("E14", "E16", "E18", "E20"))
plot_bars(tshifted_data,
          title = "Fractional abundance of glycans in temperature shifted experiments")

ggsave(filename = "figures/corrected_frac_ab_barplot_tshifted_experiments_4tp.png",
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)
# # prepare data for differential analysis ----------------------------------
# 
# data.matrix <- corr_abundance_data %>%
#   select(glycoform1, corr_abundance, condition_br_tp) %>%
#   pivot_wider(values_from = corr_abundance,
#               names_from = condition_br_tp) %>%
#   column_to_rownames('glycoform1') %>%
#   as.matrix()
# 
# meta <- tibble(sample_name = colnames(data.matrix)) %>%
#   separate(col = sample_name,
#            into = c('condition', 'br', 'timepoint'),
#            sep = "_",
#            remove = FALSE
#           )
# 
# save(data.matrix, meta, file = "analysis/matrix_meta.RData")
# 
# # prepare data for mirror plots -------------------------------------------------
# 
# make_wider_table <- function(data,
#                              condition = "[ABC]"
# ){
#   tp_264 <- data %>%
#     filter(grepl(condition, condition_br_tp)) %>%
#     filter(grepl("264", condition_br_tp)) %>%
#     rename(frac_abundance_264 = corr_abundance,error_264 = corr_abundance_error)
#   
#   tp_120 <- data %>%
#     filter(grepl(condition, condition_br_tp)) %>%
#     filter(grepl("120", condition_br_tp)) %>%
#     rename(frac_abundance_120 = corr_abundance,error_120 = corr_abundance_error) %>% 
#     select(frac_abundance_120,error_120) %>%
#     mutate(frac_abundance_120 = -frac_abundance_120)
#   
#   table_wider <- tp_264 %>%
#     cbind(tp_120) %>%
#     mutate(condition_br = str_extract(condition_br_tp, "([^_]+_[^_]+)")) %>%
#     mutate(condition = str_extract(condition_br_tp, "([^_]+)")) %>%
#     rename(modcom_name = glycoform1) %>%
#     mutate(condition = factor(condition, levels = c("F","E","D","C","B","A")))
#   
#   # return(table_wider)
#   table_wider
#   
# }
# 
# abc_wider <- make_wider_table(data = corr_abundance_data, 
#                               condition = "[ABC]")
# ab_wider <- make_wider_table(data = corr_abundance_data, 
#                               condition = "[AB]")
# ad_wider <- make_wider_table(data = corr_abundance_data, 
#                              condition = "[AD]")
# bc_wider <- make_wider_table(data = corr_abundance_data, 
#                              condition = "[BC]")
# de_wider <- make_wider_table(data = corr_abundance_data, 
#                              condition = "[DE]")
# ef_wider <- make_wider_table(data = corr_abundance_data, 
#                              condition = "[EF]")
# def_wider <- make_wider_table(data = corr_abundance_data, 
#                               condition = "[DEF]")
# abcdef_wider <- make_wider_table(data = corr_abundance_data, 
#                               condition = "[ABCDEF]")
# 
# #plot mirror plot --------------------------------------------------------------
# plot_mirror_plot <- function(data,
#                              condition = "ABC",
#                              grouping_var = condition_br,
#                              legend_cols = 3,
#                              color_mapping = color_mapping,
#                              title = "plot_title") {
#   
#   # # Define y-axis breaks and labels
#   y_breaks <- c(-60, -30, 0, 30, 60)
# 
#   ggplot(data, aes(x = modcom_name)) +
#     geom_col(aes(y = frac_abundance_264, fill = {{grouping_var}}), 
#              position = position_dodge(width = 0.9)) +
#     geom_col(aes(y = frac_abundance_120, fill = {{grouping_var}}), 
#              position = position_dodge(width = 0.9)) +
#     geom_errorbar(
#       aes(
#         ymin = frac_abundance_264 - error_264,
#         ymax = frac_abundance_264 + error_264,
#         group = {{grouping_var}}
#       ),
#       position = position_dodge(.9),
#       width = .5,
#       linewidth = .25
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = frac_abundance_120 - error_120,
#         ymax = frac_abundance_120 + error_120,
#         group = {{grouping_var}}
#       ),
#       position = position_dodge(.9),
#       width = .5,
#       linewidth = .25
#     ) +
#     geom_smooth(
#       aes(y = frac_abundance_264, group = {{grouping_var}}, color = {{grouping_var}}),
#       method = "loess", span = 0.5, se = FALSE
#     ) +
#     geom_smooth(
#       aes(y = frac_abundance_120, group = {{grouping_var}}, color = {{grouping_var}}),
#       method = "loess", span = 0.5, se = FALSE
#     ) +
#     scale_fill_manual(values = color_mapping, 
#                       breaks = names(color_mapping)) +
#     scale_color_manual(values = color_mapping, 
#                        breaks = names(color_mapping)) +
#     scale_y_continuous(name = "fractional abundance (%)",
#                        breaks = y_breaks, 
#                        labels = \(x) abs(x), 
#                        limits = c(-60,60)) +
#     coord_flip() +
#     xlab("") +
#     geom_hline(yintercept = 0, linewidth = .35) +
#     guides(fill = guide_legend(ncol = legend_cols)) +
#     annotate("text", x = 10, y = -50, label = "120 hours", alpha = 0.5) +
#     annotate("text", x = 10, y = 50, label = "264 hours", alpha = 0.5) +
#     theme_bw() +
#     theme(text = element_text(size = 10, 
#                               face = "plain",
#                               family = "sans"),
#           axis.text.y = element_text(colour = "black", hjust = 0.5),
#           axis.text = element_text(colour = "black"),
#           axis.ticks.y = element_blank(),
#           legend.title = element_blank(),
#           legend.text = element_text(size = 10),
#           legend.key.height = unit(0.3, 'cm'),
#           legend.key.width = unit(0.3, 'cm'),
#           legend.position = "bottom",
#           panel.border = element_blank(),
#           panel.grid.major.y = element_blank(),
#           panel.grid.minor = element_blank(),
#           ) +
#     ggtitle(title) +
#   NULL  
# 
#   ggsave(filename = paste0("figures/corrected_frac_ab_mirror_barplot_",condition,"_LH_colorscale.png"),
#          height = 180,
#          width = 180,
#          units = "mm",
#          dpi = 600)
# }
# 
# # Define the colors
# color_mapping_abcdef_br <- c(
#   "A_2" = "#EE3377",
#   "A_3" = "#EE3377",
#   "B_2" = "#56B4E9",
#   "C_1" = "#009E73",
#   "C_2" = "#009E73",
#   "C_3" = "#009E73",
#   "D_1" = "#CC79A7",
#   "D_2" = "#CC79A7",
#   "D_3" = "#CC79A7",
#   "E_1" = "#EE7631",
#   "E_2" = "#EE7631",
#   "E_3" = "#EE7631",
#   "F_1" = "#0072B2",
#   "F_2" = "#0072B2",
#   "F_3" = "#0072B2"
# )
# 
# 
# #currently not working, because color is mapped to condition only, while grouping is by condition_br
# plot_mirror_plot(data = abc_wider,
#                  condition = "ABC",
#                  grouping_var = abc_wider$condition_br,
#                  color_mapping = color_mapping_abcdef_br)
# 
# plot_mirror_plot(data = def_wider,
#                  condition = "DEF",
#                  grouping_var = def_wider$condition_br,
#                  color_mapping = color_mapping_abcdef_br)
# 
# plot_mirror_plot(data = ab_wider,
#                  condition = "AB",
#                  grouping_var = ab_wider$condition_br,
#                  color_mapping = color_mapping_abcdef_br,
#                  title = "Impact of Gal+")
# 
# plot_mirror_plot(data = ad_wider,
#                  condition = "AD",
#                  grouping_var = ad_wider$condition_br,
#                  color_mapping = color_mapping_abcdef_br,
#                  title = "24h_vs_48h feeding interval")
# 
# plot_mirror_plot(data = bc_wider,
#                  condition = "BC",
#                  grouping_var = bc_wider$condition_br,
#                  color_mapping = color_mapping_abcdef_br,
#                  title = "8g_vs_4g glucose every 48h")
# 
# plot_mirror_plot(data = de_wider,
#                  condition = "DE",
#                  grouping_var = de_wider$condition_br,
#                  color_mapping = color_mapping_abcdef_br,
#                  title = "HIPDOG impact")
# 
# plot_mirror_plot(data = ef_wider,
#                  condition = "EF",
#                  grouping_var = ef_wider$condition_br,
#                  color_mapping = color_mapping_abcdef_br,
#                  title = "Gal+ effect on HIPDOG")
# 
# # plots conditions compare ------------------------------------------------
# 
# # Filter rows where the second part of the string contains "2"
# filtered_corr_abundance_data <- corr_abundance_data %>%
#   filter(str_detect(condition_br_tp, "^[^_]*_2_"))
# 
# filt_abcdef_wider <- make_wider_table(data = filtered_corr_abundance_data, 
#                               condition = "[ABCDEF]")
# filt_abc_wider <- make_wider_table(data = filtered_corr_abundance_data, 
#                                 condition = "[ABC]")
# filt_def_wider <- make_wider_table(data = filtered_corr_abundance_data, 
#                                 condition = "[DEF]")
# 
# # Define the colors
# color_mapping_condition <- c(
#   "A" = "#EE3377",
#   "B" = "#56B4E9",
#   "C" = "#009E73",
#   "D" = "#CC79A7",
#   "E" = "#EE7631",
#   "F" = "#0072B2"
# )
# 
# plot_mirror_plot(data = filt_abcdef_wider,
#                  condition = "ABCDEF_cond_comp",
#                  grouping_var = filt_abcdef_wider$condition,
#                  legend_cols = 2,
#                  color_mapping = color_mapping_condition)
# 
# plot_mirror_plot(data = filt_abc_wider,
#                  condition = "ABC_cond_comp",
#                  grouping_var = condition,
#                  legend_cols = 3,
#                  color_mapping = color_mapping_condition)
# 
# plot_mirror_plot(data = filt_def_wider,
#                  condition = "DEF_cond_comp",
#                  grouping_var = filt_def_wider$condition,
#                  legend_cols = 3,
#                  color_mapping = color_mapping_condition)
# 
# 
# # plot some heatmaps ------------------------------------------------------
# # plot data as a heatmap --------------------------------------------------
# glycan_order_120 <- c(6,7, 8, 5, 9, 4, 1, 2, 10, 3)
# glycan_order_264 <- c(6, 7, 5, 8, 1, 4, 9, 2, 10, 3)
# 
# plot_heatmap <-  function(data,
#                           timepoint = "120",
#                           glycan_order = c(6, 7, 8, 5, 9, 4, 1, 2, 10, 3)) {
#   
#   data.matrix <- data %>%
#     select(glycoform1, corr_abundance, condition_br_tp) %>%
#     filter(grepl(timepoint, condition_br_tp)) %>%
#     pivot_wider(names_from = condition_br_tp, values_from = corr_abundance) %>%
#     rename_with(~ str_replace(., paste0("_",timepoint), "")) %>%
#     column_to_rownames("glycoform1") %>%
#     arrange(glycan_order)  %>%
#     as.matrix() 
#   
#   scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row 
#   
#   #heatmap settings
#   BASE_TEXT_SIZE_PT <- 9
#   ht_opt(
#     simple_anno_size = unit(1.5, "mm"),
#     COLUMN_ANNO_PADDING = unit(1, "pt"),
#     DENDROGRAM_PADDING = unit(1, "pt"),
#     HEATMAP_LEGEND_PADDING = unit(1, "mm"),
#     ROW_ANNO_PADDING = unit(1, "pt"),
#     TITLE_PADDING = unit(2, "mm"),
#     heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
#     heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
#     heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
#     heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
#     legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
#     legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
#     legend_border = FALSE
#   )
#   #color scheme
#   f1 = colorRamp2(seq(-max(abs(scaled.data.matrix)),
#                       max(abs(scaled.data.matrix)),
#                       length = 9),
#                   c("seagreen4",
#                     "seagreen3",
#                     "seagreen2",
#                     "seagreen1",
#                     "gold",
#                     "darkorchid1",
#                     "darkorchid2",
#                     "darkorchid3",
#                     "darkorchid4"),
#                   space = "RGB")
#   
#   png(filename = paste0("figures/heatmap_",timepoint,".png"),    
#       height = 9,
#       width = 9,
#       units = "cm",
#       res = 600)
#   
#   draw(Heatmap(scaled.data.matrix,
#                col = f1,
#                cluster_rows = FALSE,
#                rect_gp = gpar(col = "white", lwd = 2),
#                name = "z-score of fractional abundance",
#                row_gap = unit(4, "pt"),
#                column_gap = unit(4, "pt"),
#                width = unit(4, "mm") * ncol(scaled.data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
#                height = unit(4, "mm") * nrow(scaled.data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
#                show_row_names = TRUE,
#                heatmap_legend_param = list(direction = "horizontal")
#   ),
#   heatmap_legend_side = "bottom")
#   
#   dev.off()
# }
# 
# data.matrix_120_func <- plot_heatmap(data = corr_abundance_data) 
# 
# data.matrix_264_func <- plot_heatmap(data = corr_abundance_data,
#                                      timepoint = "264",
#                                      glycan_order = glycan_order_264) 
# 
# 
# # colorscheme for conditions ----------------------------------------------
# 
# # Define the colors
# color_mapping <- c(
#   "A" = "#EE3377",
#   "B" = "#56B4E9",
#   "C" = "#009E73",
#   "D" = "#CC79A7",
#   "E" = "#EE7631",
#   "F" = "#0072B2"
# )
# 
# 
# # plot glycans individually, bar plots ------------------------------------
# 
# filtered_tp_264 <- filtered_corr_abundance_data %>%
#   filter(grepl("264", condition_br_tp)) %>%
#   mutate(condition = str_extract(condition_br_tp, "([^_]+)")) 
#  
# filtered_tp_120 <- filtered_corr_abundance_data %>%
#   filter(grepl("120", condition_br_tp)) %>%
#   mutate(condition = str_extract(condition_br_tp, "([^_]+)")) 
# 
# 
# plot_faceted <- function(data_to_plot){
#   ggplot(data_to_plot, aes(x = condition)) +
#     geom_col(aes(y = corr_abundance, fill = condition), 
#              position = position_dodge(width = 0.9)) +
#     geom_errorbar(
#       aes(
#         ymin = corr_abundance - corr_abundance_error,
#         ymax = corr_abundance + corr_abundance_error,
#         group = condition
#       ),
#       position = position_dodge(.9),
#       width = .5,
#       linewidth = .25
#     ) + 
#     facet_wrap(~glycoform1) +
#     scale_fill_manual(values = color_mapping, 
#                       breaks = names(color_mapping)) +
#     xlab("")
# 
# }
# 
# plot_faceted(data_to_plot = filtered_tp_120)  
# plot_faceted(data_to_plot = filtered_tp_264) 
