library(tidyverse)
library(janitor)
# percent of proteins with missed cleavages -------------------------------


# protein coverage --------------------------------------------------------
subfolders <- c(
              "ptm_E13_E15_trypsin.presets",
              "ptm_E16_E18_trypsin.presets",
              "ptm_E19_E20_TB_KB_trypsin.presets"
               )

# Initialize an empty list to store data frames
protein_coverage_list <- list()

for (subfolder in subfolders) {
  data_path <- paste0("analysis/peptide_mapping/", 
                    subfolder, 
                    "/protein_coverage_",
                    subfolder,
                    ".csv")

  protein_coverage <-  read_csv(data_path) %>%
    clean_names() %>%
    filter(!row_number() %in% c(1, 2)) %>%
    separate(sample_name,
            into = c("date", 
                      "initials",
                      "id1",
                      "id2",
                      "column_details",
                      "experiment",
                      "timepoint",
                      "digest",
                      "technical_replicate",
                      "acquisition_number"),
            sep = "_",
            remove = FALSE)

  stats_protein_coverage <- protein_coverage %>%
    group_by(experiment, timepoint, protein_name) %>%
    summarise(mean_percent = mean(protein_coverage_percent),
              sd_percent = sd(protein_coverage_percent)) %>%
    mutate(timepoint = ifelse(timepoint == "500ng" | timepoint == "T", "harvest", timepoint),
          experiment = ifelse(experiment == "NISTMAb" , "NISTMAB_TB", experiment),
          experiment = ifelse(experiment == "Trypsin" , "NISTMAB_KB", experiment),
          experiment_timepoint = paste(experiment, timepoint, sep = "_"))

  protein_coverage_list[[subfolder]] <- stats_protein_coverage
}

# Combine all data frames into a single one
protein_coverage_all <- bind_rows(protein_coverage_list, .id = "subfolder")


write_csv(x = protein_coverage_all,
          file = "analysis/peptide_mapping/protein_coverage.csv")

# plot the data -----------------------------------------------------------
color_mapping_experiment <- c(
  "E13" = "#FD8D3C",
  "E14" = "#9E9AC8",
  "E15" = "#F16913",
  "E16" = "#807DBA",
  "E17" = "#D94801",
  "E18" = "#6A51A3",
  "E19" = "#A63603",
  "E20" = "#54278F",
  "NISTMAB_TB" = "#1b9e77",
  "NISTMAB_KB" = "#e7298a"
)

ggplot(protein_coverage_all, aes(x = experiment_timepoint, y = mean_percent, fill = experiment)) +
  geom_col(
    position = position_dodge(width = 0.9)
  ) + 
  geom_hline(yintercept = 60, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_errorbar(
    aes(
      ymin = mean_percent - sd_percent,
      ymax = mean_percent + sd_percent,
      group = experiment_timepoint
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_manual(values = color_mapping_experiment) +
  scale_y_continuous(name = "protein_coverage [%]",
                     limits = c(0,100)) +
  facet_wrap(~protein_name, ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  )

ggsave("figures/peptide_mapping/protein_coverage.png",
       width = 150,
       height = 150,
       units = "mm",
       dpi = 300)


data_to_plot <- protein_coverage_all %>% filter((experiment %in% c("NISTMAB_TB", "NISTMAB_KB")))

ggplot(data_to_plot, aes(x = experiment_timepoint, y = mean_percent, fill = experiment)) +
  geom_col(
    position = position_dodge(width = 0.9)
  ) + 
  geom_hline(yintercept = 60, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_errorbar(
    aes(
      ymin = mean_percent - sd_percent,
      ymax = mean_percent + sd_percent,
      group = experiment_timepoint
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_manual(values = color_mapping_experiment) +
  scale_y_continuous(name = "protein_coverage [%]",
                     limits = c(0,100)) +
  scale_x_discrete(name = "") +
  facet_wrap(~protein_name, ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "none"
  )

ggsave("figures/peptide_mapping/protein_coverage_reference.png",
       width = 50,
       height = 150,
       units = "mm",
       dpi = 300)
