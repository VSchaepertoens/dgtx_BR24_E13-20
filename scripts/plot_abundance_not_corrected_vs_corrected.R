library(tidyverse)

load("analysis/corr_abundance_data.RData")
load("analysis/500_ppm/abundance_data_none.RData")

abundance_data_averaged <- abundance_data_averaged %>% 
  mutate(hexose_bias = "not_corrected") %>%
  mutate(modcom_name = str_replace_all(modcom_name, "A2", ""))  

corr_abundance_data <- corr_abundance_data %>% 
  mutate(hexose_bias = "corrected") %>%
  rename(modcom_name = glycoform1, 
         frac_abundance = corr_abundance, 
         error = corr_abundance_error) %>%
  mutate(modcom_name = str_replace_all(modcom_name, c("G0F/G1F" = "G1F/G0F",
                                                      "G1F/G2F" = "G2F/G1F")))

merged_data <- abundance_data_averaged %>%
  rbind(corr_abundance_data) %>%
  mutate(experiment_hb = paste(experiment, hexose_bias, sep = "_")) %>%
  mutate(modcom_name = factor(modcom_name, levels = c("none/G0F",
                                                    "none/G1F",
                                                    "none/G2F",
                                                    "G0/G0",
                                                    "G0/G0F",
                                                    "G0F/G0F",
                                                    "G1F/G0F",
                                                    "G1F/G1F",
                                                    "G2F/G1F",
                                                    "G2F/G2F")
  ))


# plot corrected vs not corrected data ------------------------------------

color_mapping_experiment <- c(
  "E13_corrected" = "#FD8D3C",       # Opaque
  "E13_not_corrected" = "#FD8D3C80", # Semi-transparent using RGBA (80 is ~50% opacity)
  "E14_corrected" = "#9E9AC8",       # Opaque
  "E14_not_corrected" = "#9E9AC880", # Semi-transparent
  "E15_corrected" = "#F16913",
  "E15_not_corrected" = "#F1691380", # Semi-transparent
  "E16_corrected" = "#807DBA",
  "E16_not_corrected" = "#807DBA80", # Semi-transparent
  "E17_corrected" = "#D94801",
  "E17_not_corrected" = "#D9480180", # Semi-transparent
  "E18_corrected" = "#6A51A3",
  "E18_not_corrected" = "#6A51A380", # Semi-transparent
  "E19_corrected" = "#A63603",
  "E19_not_corrected" = "#A6360380", # Semi-transparent
  "E20_corrected" = "#54278F",
  "E20_not_corrected" = "#54278F80"  # Semi-transparent
)

plot_bars <- function(data,
                      title = "Fractional abundance",
                      row_number = 1,
                      legend_row_number = 1) {
  
  ggplot(data, aes(x = modcom_name, y = frac_abundance)) +
    
    # Main bar plot with fill color based on experiment and alpha for uncorrected
    geom_col(aes(
      fill = experiment_hb,  # Color based on experiment
      alpha = ifelse(hexose_bias == "corrected", 1, 0.5)  # Alpha for uncorrected
    ), position = position_dodge(width = 0.9)) +
    
    # Ensure alpha is handled correctly
    scale_alpha_identity() +
    
    # Error bars with dodging
    geom_errorbar(
      aes(
        ymin = frac_abundance - error,
        ymax = frac_abundance + error,
        group = experiment_hb
      ),
      position = position_dodge(.9),
      width = .5,
      linewidth = .25
    ) +
    
    # Customize legend
    guides(fill = guide_legend(nrow = legend_row_number)) +
    
    # Facet by timepoint
    facet_wrap(~ timepoint, nrow = row_number) +
    
    # Use the color mapping based on experiments
    scale_fill_manual(values = color_mapping_experiment) +
    
    # Y-axis limits and label
    scale_y_continuous(name = "Fractional abundance (%)",
                       limits = c(0, 65)) +
    
    # X-axis and general theme settings
    xlab("") +
    theme_bw() +
    theme(
      text = element_text(size = 10, face = "bold", family = "sans"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text = element_text(colour = "black"),
      # legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.border = element_blank()
    ) +
    
    # Title
    ggtitle(title) +
    NULL
}


plot_bars(merged_data %>%
            filter(experiment %in% c("E20") & timepoint %in% c("120", "216", "288", "336")), 
          row_number = 1,
          legend_row_number = 2)

ggsave(filename = "figures/cafog_comparison_e20_6times30.png",
       height = 6,
       width = 30,
       units = "cm",
       dpi = 600)
