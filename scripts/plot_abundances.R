library(RColorBrewer)
library(tidyverse)


# define analysis of pngase F digested or not digested data ---------------

pngase <- "none" # "none"

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv(paste0("analysis/overview_",pngase,"_merged.csv")) %>%
  # filter(filename != "20240904_TB_cnistcho_E13_192_none_4_29.mzML") %>%#possibly an outlier
  {.}

# load abundances using a for loop  ---------------------------------------

abundance_data <- NULL
for (i in 1:nrow(samples_table)) {
  file_path <- paste0(samples_table[i, "analysis_path"], "/frac_ab_tb_cs50.csv")
  
  abundance_data <- rbind(abundance_data,
                          read_csv(file_path)
                          )
}

abundance_data <- abundance_data %>%
  separate(file_name,
           into = c("data_ymd", 
                    "initials",
                    "product",
                    "experiment",
                    # "biological_replicate",
                    "timepoint", 
                    "pngase",
                    "technical_replicate",
                    "acquisition_number"
                    ),
           sep = "_",
           remove = FALSE) %>%
  mutate(experiment_tp = paste(experiment,timepoint, sep = "_"))
  

# calculate mean and sd and plot ------------------------------------------
if (pngase == "none") {
  modcom_levels <- c("none/A2G0F",
                     "none/A2G1F",
                     "none/A2G2F",
                     "A2G0/A2G0",
                     "A2G0/A2G0F",
                     "A2G0F/A2G0F",
                     "A2G1F/A2G0F",
                     "A2G1F/A2G1F",
                     "A2G2F/A2G1F",
                     "A2G2F/A2G2F")
} else if (pngase == "pngase")  {
  modcom_levels <- c("3xHex", "2xHex","1xHex","none")
}

abundance_data_averaged <- abundance_data %>% 
  group_by(modcom_name, experiment_tp) %>%
  summarise(frac_abundance = mean(frac_ab),
            error = sd(frac_ab)) %>%
  mutate(modcom_name = factor(modcom_name, levels = modcom_levels)) %>%
  ungroup() %>%
  separate(experiment_tp,
           into = c("experiment",
                    "timepoint"
                    ),
           sep = "_",
           remove = FALSE)

save(abundance_data,
     abundance_data_averaged,
     file = paste0("analysis/abundance_data_",pngase,".RData"))

# plot char runs data -----------------------------------------------------

ggplot(abundance_data_averaged, aes(x = modcom_name, y = frac_abundance, fill = experiment)) +
  geom_col(
    position = position_dodge(width = 0.9)  
  ) + 
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
      group = experiment_tp
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  facet_wrap(~ timepoint) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))


# # vertical bar plot -------------------------------------------------------
# plot_vertical_barplot <- function(data_to_plot,
#                                   condition = "A",
#                                   timepoint = "120"
#                                   ){
# 
#   data_to_plot <- data_to_plot %>%
#     filter(grepl(condition, condition_br_tp)) %>%
#     filter(grepl(timepoint, condition_br_tp))
#     
#   data_to_plot %>%
#     ggplot(aes(x = modcom_name, y = frac_abundance, fill = condition_br_tp)) +
#     geom_col(
#       position = position_dodge(width = 0.9)  
#     ) +
#     geom_errorbar(
#       aes(
#         ymin = frac_abundance - error,
#         ymax = frac_abundance + error,
#         group = condition_br_tp
#       ),
#       position = position_dodge(.9),
#       width = .5,
#       linewidth = .25
#     ) +
#     xlab("") +
#     ylim(0, 100) +
#     ylab("fractional abundance (%)") +
#     geom_hline(yintercept = 0, linewidth = .35) +
#     coord_flip() +
#     theme_bw() +
#     guides(fill = guide_legend(ncol = 3)) +
#     theme(text = element_text(size = 9, 
#                               # face = "bold", 
#                               family = "sans"),
#           axis.text.y = element_text(colour = "black", hjust = 0.5),
#           axis.text = element_text(colour = "black"),
#           axis.ticks.y = element_blank(),
#           legend.title = element_blank(),
#           legend.text = element_text(size = 9),
#           legend.key.height = unit(0.3, 'cm'),
#           legend.key.width = unit(0.3, 'cm'),
#           legend.position = "bottom",
#           panel.border = element_blank(),
#           panel.grid.major.y = element_blank(),
#           panel.grid.minor = element_blank(),
#     ) 
#   ggsave(filename = paste0("figures/frac_ab_barplot",condition,"_",timepoint,"_",pngase,".png"),    
#          height = 160,
#          width = 160,
#          units = "mm",
#          dpi = 600)
# }
# 
# plot_vertical_barplot(abundance_data_averaged, 
#                       condition = "[ABC]",
#                       timepoint = "120") 
# 
# 
# # mirror plots --------------------------------------------------------------
# 
# make_wider_table <- function(data,
#                              condition = "[ABC]"
#                              ){
# tp_120 <- data %>%
#   filter(grepl(condition, condition_br_tp)) %>%
#   filter(grepl("120", condition_br_tp)) %>%
#   rename(frac_abundance_120 = frac_abundance,error_120 = error)
# 
# tp_264 <- data %>%
#   filter(grepl(condition, condition_br_tp)) %>%
#   filter(grepl("264", condition_br_tp)) %>%
#   rename(frac_abundance_264 = frac_abundance,error_264 = error) %>% 
#   select(frac_abundance_264,error_264) %>%
#   mutate(frac_abundance_264 = -frac_abundance_264)
# 
# table_wider <- tp_120 %>%
#   cbind(tp_264) %>%
#   mutate(condition_br = str_extract(condition_br_tp, "([^_]+_[^_]+)"))
# 
# return(table_wider)
# }
# 
# abc_wider <- make_wider_table(data = abundance_data_averaged, condition = "[ABC]")
# def_wider <- make_wider_table(data = abundance_data_averaged, condition = "[DEF]")
# 
# #plot mirror plot
# ggplot(def_wider, aes(x = modcom_name)) +
#   geom_col(aes(y = frac_abundance_264, fill = condition_br), position = position_dodge(width = 0.9)) +
#   geom_col(aes(y = frac_abundance_120, fill = condition_br), position = position_dodge(width = 0.9)) +
#   geom_errorbar(
#     aes(
#       ymin = frac_abundance_264 - error_264,
#       ymax = frac_abundance_264 + error_264,
#       group = condition_br
#     ),
#     position = position_dodge(.9),
#     width = .5,
#     linewidth = .25
#   ) +
#   geom_errorbar(
#     aes(
#       ymin = frac_abundance_120 - error_120,
#       ymax = frac_abundance_120 + error_120,
#       group = condition_br_tp
#     ),
#     position = position_dodge(.9),
#     width = .5,
#     linewidth = .25
#   ) +
#   coord_flip() +
#   ylim(-60, 60) +
#   xlab("") +
#   ylab("fractional abundance (%)") +
#   geom_hline(yintercept = 0, linewidth = .35) +
#   coord_flip() +
#   theme_bw() +
#   guides(fill = guide_legend(ncol = 3)) +
#   theme(text = element_text(size = 9, 
#                             # face = "bold", 
#                             family = "sans"),
#         axis.text.y = element_text(colour = "black", hjust = 0.5),
#         axis.text = element_text(colour = "black"),
#         axis.ticks.y = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 9),
#         legend.key.height = unit(0.3, 'cm'),
#         legend.key.width = unit(0.3, 'cm'),
#         legend.position = "bottom",
#         panel.border = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#   ) 
# 
# ggsave(filename = paste0("figures/frac_ab_mirror_barplot_ABC.png"),    
#        height = 160,
#        width = 160,
#        units = "mm",
#        dpi = 600)
