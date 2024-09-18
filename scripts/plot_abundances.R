library(RColorBrewer)
library(tidyverse)


# define analysis of pngase F digested or not digested data ---------------

pngase <- "pngase" # "none" or "pngase"

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
           remove = FALSE) %>%
  mutate(experiment = factor(experiment, 
                           levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20")
                           )) %>%
  filter(experiment != c("E17")) %>%
  filter(timepoint %in% c("120", "216", "288","336")) %>%
  {.}

save(abundance_data,
     abundance_data_averaged,
     file = paste0("analysis/abundance_data_",pngase,".RData"))

load(paste0("analysis/abundance_data_",pngase,".RData"))
# plot char runs data -----------------------------------------------------
# Define the colors
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


plot_bars <- function(data){
ggplot(data, aes(x = modcom_name, y = frac_abundance, fill = experiment)) +
  geom_col(
    position = position_dodge(width = 0.9)  
  ) + 
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
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
                       # limits = c(0,65)
                       ) +
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
    NULL
}

plot_bars(abundance_data_averaged)

ggsave(filename = "figures/pngase_frac_ab_barplot_minus17_all_experiments_4tp.png",
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)

constant_data <- abundance_data_averaged %>%
  filter(experiment %in% c("E13", "E15", "E17", "E19"))
plot_bars(constant_data)

ggsave(filename = "figures/frac_ab_barplot_constant.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600)

tshifted_data <- abundance_data_averaged %>%
  filter(experiment %in% c("E14", "E16", "E18", "E20"))
plot_bars(tshifted_data)

ggsave(filename = "figures/frac_ab_barplot_tshifted.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600)

e13e15_data <- abundance_data_averaged %>%
  filter(experiment %in% c("E13", "E15"))
plot_bars(e13e15_data)

ggsave(filename = "figures/frac_ab_barplot_e13e15.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600)

e18e20_data <- abundance_data_averaged %>%
  filter(experiment %in% c("E18", "E20"))
plot_bars(e18e20_data)

ggsave(filename = "figures/frac_ab_barplot_e18_e20.png",
         height = 100,
         width = 200,
         units = "mm",
         dpi = 600)

e15e16_data <- abundance_data_averaged %>%
  filter(experiment %in% c("E15", "E16"))
plot_bars(e15e16_data)

ggsave(filename = "figures/frac_ab_barplot_e15_e16.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600)

e13e14_data <- abundance_data_averaged %>%
  filter(experiment %in% c("E13", "E14"))
plot_bars(e13e14_data)

ggsave(filename = "figures/frac_ab_barplot_e13_e14.png",
       height = 100,
       width = 200,
       units = "mm",
       dpi = 600)

# plot line plots for each glycan seprately -------------------------------
# Check for missing values
sum(is.na(abundance_data_averaged$timepoint))
sum(is.na(abundance_data_averaged$frac_abundance))
sum(is.na(abundance_data_averaged$experiment))

str(abundance_data_averaged) # Check the structure of your dataframe
abundance_data_averaged$timepoint <- as.numeric(as.character(abundance_data_averaged$timepoint))


ggplot(abundance_data_averaged, aes(x = timepoint, y = frac_abundance, color = experiment)) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) + # Remove fullrange = TRUE
  facet_wrap(~modcom_name, 
             scales = "free_y",
             nrow = 2)

ggplot(abundance_data_averaged %>% filter(experiment %in% c("E14")), 
       aes(x = timepoint, y = frac_abundance, fill = experiment)) +
  geom_col(
    position = position_dodge(width = 0.9)  
  ) + 
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
      group = experiment
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_manual(values = color_mapping_experiment, 
                    breaks = names(color_mapping_experiment)) +
  facet_wrap(~modcom_name, 
             scales = "free_y",
             nrow = 2)


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
