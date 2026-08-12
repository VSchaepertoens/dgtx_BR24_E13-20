library(RColorBrewer)
library(tidyverse)


# define analysis of pngase F digested or not digested data ---------------

pngase <- "none" # "none" or "pngase"

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv(paste0("analysis/overview_",pngase,"_merged.csv")) %>%
  filter(filename.x != "20251002_TB_cnistcho_E13_312_none_1_336.mzML") %>% #possibly an outlier
  filter(filename.x != "20251113_TB_cnistcho_E20_288_none_3_1368.mzML") %>% #possibly an outlier
  # filter(filename.x != "20251113_TB_cnistcho_E20_288_none_2_1367.mzML") %>% #possibly an outlier
  
  # filter(sample_name != "20251002_TB_Nistmab_150mg_l_pngase") %>% # nistmab control
  
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
                           # levels = c("E15", "E17", "E19", "E16", "E18", "E20")
                           ),
         timepoint = factor(timepoint, 
                             # levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20")
                             levels = c("72", "96", "120", "144", "168", "192","216",
                                        "240", "264", "288", "312", "336", "360"))
         ) %>%
  {.}

save(abundance_data,
     abundance_data_averaged,
     file = paste0("analysis/abundance_data_",pngase,".RData"))

load(paste0("analysis/abundance_data_",pngase,".RData"))
# plot char runs data -----------------------------------------------------
# Define the colors
color_mapping_experiment <- c(
  "E13" = "#FD8D3C",
  "E14" = "#9E9AC8",
  "E15" = "#F16913",
  "E16" = "#807DBA",
  "E17" = "#D94801",
  "E18" = "#6A51A3",
  "E19" = "#A63603",
  "E20" = "#54278F"
)


plot_bars <- function(data,
                      row_number = 1){
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
    guides(fill = guide_legend(nrow = row_number)) +
    facet_wrap(~ timepoint, nrow = row_number) +
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

abundance_data %>%
  filter(experiment_tp %in% "E20_288")

plot_bars(abundance_data_averaged, row_number = 2)

ggsave(filename = "figures/20251001_TB_cNISTCHO_CharRUns_rem_allTP/pngase_frac_ab_barplot_all_experiments.png",
       height = 200,
       width = 250,
       units = "mm",
       dpi = 600)

plot_bars(abundance_data_averaged %>% filter(experiment != c("E17")), 
          row_number = 2)

ggsave(filename = "figures/pngase_frac_ab_barplot_minus17_all_experiments.png",
       height = 200,
       width = 250,
       units = "mm",
       dpi = 600)

plot_bars(abundance_data_averaged %>% filter(timepoint %in% c("120", "216", "288","336") & experiment != c("E17")), 
          row_number = 1)

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


