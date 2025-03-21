# script to load all peptide mapping data and merge across samples and replicates
library(tidyverse)
library(janitor)
library(circlize)
library(ComplexHeatmap)
library(ggpattern)

# load data -------------------------------------------------------
subfolders <- c(
                "E17_originator_TB_KB_trypsin",
                "ptm_E13toE18_originator_TB_KB_trypsin.presets",
                "ptm_E19_E20"
                )

# stats for individual sample -------------------------------
for (subfolder in subfolders) {
  data_path <- paste0("analysis/peptide_mapping/", 
                      subfolder, 
                      "/nglycans_results.csv")
  print(data_path)

  data <- read_csv(data_path) %>%
    select(
      "MS\r\nAlias name",
      "Group",
      "Protein\r\nname",
      "Sequence",
      "Mod.\r\nSummary",
      "Glycans",
      "Labels",
      "z",
      "Total XIC AUC\r\nAveragine",
      "Validate",
      "Digest\r\nname"
    ) %>%
    clean_names() %>%
    mutate(glycans = case_when(
      labels == "unglycosylated" ~ "unglycosylated",
      TRUE ~ glycans
    ))

# sum up xic averagine for multiple charge states
data_summed <- data %>%
  group_by(ms_alias_name, glycans, labels) %>%
  summarize(sum_xic_averagine = sum(total_xic_auc_averagine)) %>%
  ungroup() %>%
  group_by(ms_alias_name) %>%
  mutate(frac_abud = sum_xic_averagine / sum(sum_xic_averagine) * 100)


#test that for a given esample, all glycan abundances sum to 100%
test <- data_summed %>%
  group_by(ms_alias_name) %>%
  summarize(frac_abud_sum = sum(frac_abud))

#expand the files names per experiment and per timepoint
data_summed <- data_summed %>%
  separate(ms_alias_name,
         into = c("date",
                  "initials1",
                  "instrument_ID",
                  "instrument_ID2",
                  "column",
                  "experiment",
                  "timepoint",
                  "digestion_enzyme",
                  "technical_replicate",
                  "acquisition_number"
                  ),
         sep = "_",
         remove = FALSE) %>%
  mutate(experiment_timepoint = paste(experiment,
                                      timepoint,
                                      sep = "_")) %>%
  {.}

# plot the abundance for a given experiment_timepoint, all replicates
p <- ggplot(data_summed, aes(x = labels, y = frac_abud)) +
  geom_col(aes(
    fill = ms_alias_name
    ), position = position_dodge(width = 0.9),
    color = "black",
    linewidth = 0.1) +
  facet_wrap(~ experiment_timepoint, nrow = 2) +
  scale_y_continuous(name = "Fractional abundance (%)",
                     limits = c(0, 80)) +
  scale_x_discrete(name = "N-glycosylation type",
                   guide = guide_axis(angle = 90)) +
  guides(fill = guide_legend(ncol = 4)) +
  ggtitle(paste0("R.EEQYnSTYR.V peptide for sample ",sub("\\..*", "", subfolder))) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 6))
plot(p)

output_data_path <- paste0("analysis/peptide_mapping/",
                           subfolder,
                           "/results_table_nglycans_REEQYnSTYRV.csv")

write_csv(data_summed,
          output_data_path)
}


# stats for biological replicates mean + sd -------------------------------

abundance_data <- NULL

for (subfolder in subfolders) {
  data_path <- paste0("analysis/peptide_mapping/", 
                      subfolder, 
                      "/results_table_nglycans_REEQYnSTYRV.csv")
  
  print(data_path)
  
  if (subfolder == "E17_originator_TB_KB_trypsin") {
    data <- read_csv(data_path) %>%
      filter(experiment == "E17")
  } else {
    data <- read_csv(data_path) 
  }
  
  abundance_data <- rbind(abundance_data,
                          data)
}

print(unique(abundance_data$experiment))

# mean + sd of biological replicates
abundance_data_summed <- abundance_data %>%
  group_by(experiment_timepoint, labels) %>%
  summarize(mean_frac_abud = mean(frac_abud),
            sd_frac_abud = sd(frac_abud)) %>%
  ungroup() %>%
  separate(experiment_timepoint,
           into = c("experiment", 
                    "timepoint"),
           sep = "_",
           remove = FALSE) %>%
  mutate(timepoint = ifelse(timepoint == "500ng" | timepoint == "T", "288", timepoint),
         experiment = ifelse(experiment == "NISTMAb" , "NISTMAB_TB", experiment),
         experiment = ifelse(experiment == "Trypsin" , "NISTMAB_KB", experiment),
         experiment_timepoint = paste(experiment, timepoint, sep = "_"))

# Define the colors from the "Paired" palette
color_mapping <- c(
    "E13_144" = "#FD8D3C",
    "E13_288" = "#FD8D3C",
    "E14_144" = "#9E9AC8",
    "E14_288" = "#9E9AC8",
    "E15_144" = "#F16913",
    "E15_288" = "#F16913",
    "E16_144" = "#807DBA",
    "E16_288" = "#807DBA",
    "E17_144" = "#D94801",
    "E17_288" = "#D94801",
    "E18_144" = "#6A51A3",
    "E18_288" = "#6A51A3",
    "E19_144" = "#A63603",
    "E19_288" = "#A63603",
    "E20_144" = "#54278F",
    "E20_288" = "#54278F",
    "NISTMAB_TB_288" = "#1b9e77",
    "NISTMAB_KB_288" = "#e7298a"
  )


# 1. per timepoint plots --------------------------------------------------

plot_per_timepoint <- function(in_data,
                               wo_originator = FALSE){
  
  in_data <- in_data %>%
    mutate(experiment_timepoint = factor(experiment_timepoint, levels = c("E13_144","E13_288", 
                                                                          "E15_144","E15_288", 
                                                                          "E17_144","E17_288",
                                                                          "E19_144","E19_288",
                                                                          "E14_144","E14_288", 
                                                                          "E16_144","E16_288", 
                                                                          "E18_144","E18_288",
                                                                          "E20_144","E20_288",
                                                                          "NISTMAB_TB_288", "NISTMAB_KB_288")))
  if (wo_originator) {
  in_data <- in_data %>% filter(!(experiment_timepoint %in% c("NISTMAB_TB_288", "NISTMAB_KB_288")))
  }
                      
  p <- ggplot(in_data, aes(x = labels, y = mean_frac_abud)) +
    geom_col(aes(fill = experiment_timepoint), 
             position = position_dodge(width = 0.9),
             color = "black",
             linewidth = 0.1) +
    geom_errorbar(
      aes(
        ymin = mean_frac_abud - sd_frac_abud,
        ymax = mean_frac_abud + sd_frac_abud,
        group = experiment_timepoint
      ),
      position = position_dodge(.9),
      width = .5,
      linewidth = .25
    ) +
    facet_wrap(~ timepoint, 
               nrow = 2,
               labeller = as_labeller(c('144' = 'timepoint 144 hours', '288' = 'timepoint 288 hours'))
               ) +
    scale_fill_manual(values = color_mapping) +
    scale_y_continuous(name = "Fractional abundance (%)") +
    scale_x_discrete(name = "N-glycosylation type") +
    ggtitle("Per timepoint N-glycans fractional abundance of peptide R.EEQYnSTYR.V") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = 'bottom',
          legend.text = element_text(size = 8))
  
  plot(p)
}

plot_per_timepoint(in_data = abundance_data_summed)
ggsave(filename = "figures/peptide_mapping/per_timepoint_EEQNSTYR_all.png",
       units = "cm",
       dpi = 300,
       width = 20,
       height = 20)
plot_per_timepoint(in_data = abundance_data_summed,
                   wo_originator = TRUE)
ggsave(filename = "figures/peptide_mapping/per_timepoint_EEQNSTYR_wo_originator.png",
       units = "cm",
       dpi = 300,
       width = 20,
       height = 20)
subset <- abundance_data_summed %>% filter(!(experiment_timepoint %in% c("E13_288", "E14_288"))) # filter out digests with many missed cleavages
plot_per_timepoint(in_data = subset)
ggsave(filename = "figures/peptide_mapping/per_timepoint_EEQNSTYR_E13_E14_288_out.png",
       units = "cm",
       dpi = 300,
       width = 25,
       height = 20)
plot_per_timepoint(in_data = subset,
                   wo_originator = TRUE)
ggsave(filename = "figures/peptide_mapping/per_timepoint_EEQNSTYR_E13_E14_288_out_wo_originator.png",
       units = "cm",
       dpi = 300,
       width = 25,
       height = 20)

# 2. per experiment plots -------------------------------------------------
plot_per_experiment <- function(in_data,
                                wo_originator = FALSE){
  
  # Filter out originator
  if (wo_originator) {
    in_data <- in_data %>% filter(!(experiment_timepoint %in% c("NISTMAB_TB_288", "NISTMAB_KB_288")))
    
  }
  
  # Precompute the pattern column
  in_data <- in_data %>%
    mutate(pattern = ifelse(str_detect(timepoint, "288"), "stripe", NA))
  
  p <- ggplot(in_data, aes(x = labels, y = mean_frac_abud)) +
    geom_col_pattern(
      aes(
        fill = experiment_timepoint,
        pattern = pattern
      ),
      position = position_dodge(width = 1),
      color = "black",
      linewidth = 0.1,
      pattern_fill = "black",  # Color of the stripes
      # pattern_density = 0.1,   # Adjust spacing of stripes
      # pattern_angle = 45       # Adjust stripe angle
    ) +
    geom_errorbar(
      aes(
        ymin = mean_frac_abud - sd_frac_abud,
        ymax = mean_frac_abud + sd_frac_abud,
        group = experiment_timepoint
      ),
      position = position_dodge(.9),
      width = .5,
      linewidth = .25
    ) +
    facet_wrap(~experiment, ncol = 2) + 
    scale_fill_manual(
      values = color_mapping
    ) +
    scale_pattern_manual(
      values = c("stripe", "none"),
      name = "Pattern Key",
      labels = c("288 Timepoint", "144 Timepoint")
    ) +
    scale_y_continuous(name = "Fractional abundance (%)") +
    scale_x_discrete(name = "N-glycosylation type") +
    guides(fill = guide_legend(ncol = 6)
           ) +
    ggtitle("Per experiment N-glycans fractional abundance of peptide R.EEQYnSTYR.V") +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    NULL
  
  plot(p)
}

plot_per_experiment(in_data = abundance_data_summed)
ggsave(filename = "figures/peptide_mapping/per_experiment_EEQNSTYR_all.png",
       units = "cm",
       dpi = 300,
       width = 20,
       height = 20)

plot_per_experiment(in_data = abundance_data_summed,
                    wo_originator = TRUE)
ggsave(filename = "figures/peptide_mapping/per_experiment_EEQNSTYR_wo_originator.png",
       units = "cm",
       dpi = 300,
       width = 20,
       height = 20)
subset <- abundance_data_summed %>% filter(!(experiment_timepoint %in% c("E13_288", "E14_288"))) # filter out digests with many missed cleavages
plot_per_experiment(in_data = subset)
ggsave(filename = "figures/peptide_mapping/per_experiment_EEQNSTYR_E13_E14_288_out.png",
       units = "cm",
       dpi = 300,
       width = 25,
       height = 20)
plot_per_experiment(in_data = subset,
                    wo_originator = TRUE)
ggsave(filename = "figures/peptide_mapping/per_experiment_EEQNSTYR_E13_E14_288_out_wo_originator.png",
       units = "cm",
       dpi = 300,
       width = 25,
       height = 20)

# 3. per glycan type plots -------------------------------------------------------------------------
plot_per_glycan <- function(in_data,
                            wo_originator = FALSE){
  
  in_data <- in_data %>%
    mutate(experiment_timepoint = factor(experiment_timepoint, levels = c("E13_144","E13_288", 
                                                                          "E15_144","E15_288", 
                                                                          "E17_144","E17_288", 
                                                                          "E19_144","E19_288",
                                                                          "E14_144","E14_288", 
                                                                          "E16_144","E16_288", 
                                                                          "E18_144","E18_288",
                                                                          "E20_144","E20_288",
                                                                          "NISTMAB_TB_288", "NISTMAB_KB_288")))
  if (wo_originator) {
    in_data <- in_data %>% filter(!(experiment_timepoint %in% c("NISTMAB_TB_288", "NISTMAB_KB_288")))
  }
  
  p <- ggplot(in_data, aes(x = experiment_timepoint, y = mean_frac_abud)) +
    geom_col(aes(fill = experiment_timepoint), 
             position = position_dodge(width = 0.9),
             color = "black",
             linewidth = 0.1) +
    geom_errorbar(
      aes(
        ymin = mean_frac_abud - sd_frac_abud,
        ymax = mean_frac_abud + sd_frac_abud,
        group = experiment_timepoint
      ),
      position = position_dodge(.9),
      width = .5,
      linewidth = .25
    ) +
    facet_wrap(~ labels, 
               ncol = 6,
               scales = "free_y"
               # labeller = as_labeller(c('144' = 'timepoint 144 hours', '288' = 'timepoint 288 hours'))
    ) +
    scale_fill_manual(values = color_mapping) +
    scale_y_continuous(name = "Fractional abundance (%)") +
    scale_x_discrete(name = "experimentNumber_timepoint") +
    ggtitle("Per N-glycan type fractional abundance of peptide R.EEQYnSTYR.V") +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  plot(p)
}

plot_per_glycan(in_data = abundance_data_summed)
ggsave(filename = "figures/peptide_mapping/per_glycan_EEQNSTYR_all.png",
       units = "cm",
       dpi = 300,
       width = 30,
       height = 20)
plot_per_glycan(in_data = abundance_data_summed, wo_originator = TRUE)
ggsave(filename = "figures/peptide_mapping/per_glycan_EEQNSTYR_wo_originator.png",
       units = "cm",
       dpi = 300,
       width = 30,
       height = 20)
subset <- abundance_data_summed %>% filter(!(experiment_timepoint %in% c("E13_288", "E14_288"))) # filter out digests with many missed cleavages
plot_per_glycan(in_data = subset)
ggsave(filename = "figures/peptide_mapping/per_glycan_EEQNSTYR_E13_E14_288_out.png",
       units = "cm",
       dpi = 300,
       width = 25,
       height = 20)
plot_per_glycan(in_data = subset,
                    wo_originator = TRUE)
ggsave(filename = "figures/peptide_mapping/per_glycans_EEQNSTYR_E13_E14_288_out_wo_originator.png",
       units = "cm",
       dpi = 300,
       width = 25,
       height = 20)

# plot data biological replicates mean as a heatmap --------------------------------------------------
data.matrix <- abundance_data_summed %>%
  select(labels, mean_frac_abud, experiment_timepoint) %>%
  pivot_wider(names_from = experiment_timepoint, values_from = mean_frac_abud) %>%
  column_to_rownames("labels") %>%
  # arrange(c(7,8, 9, 10, 11, 12, 4, 3, 6, 5, 1, 2)) %>%
  as.matrix()

## calculate z-score & plot heatmap -------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])


BASE_TEXT_SIZE_PT <- 9
ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(2, "mm"),
  heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_border = FALSE
)

#set the correct color scheme
# min(scaled.data.matrix)
# max(scaled.data.matrix)
# f1 = colorRamp2(seq(-max(abs(data.matrix), na.rm = TRUE),
#                     max(abs(data.matrix), na.rm = TRUE),
#                     length = 9),
#                 c("seagreen4",
#                   "seagreen3",
#                   "seagreen2",
#                   "seagreen1",
#                   "gold",
#                   "darkorchid1",
#                   "darkorchid2",
#                   "darkorchid3",
#                   "darkorchid4"),
#                 space = "RGB")
#set the correct color scheme
# png(filename = "figures/Jan_2024/heatmap_scaled_cafog_corrected_reordered_SC.png",    
#     height = 9,
#     width = 8.89,
#     units = "cm",
#     res = 600)


draw(Heatmap(scaled.data.matrix,
             col = rev(rainbow(10)),
             cluster_rows = FALSE,
             rect_gp = gpar(col = "white", lwd = 2),
             name = "fractional abundance",
             row_gap = unit(4, "pt"),
             column_gap = unit(4, "pt"),
             width = unit(4, "mm") * ncol(scaled.data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
             height = unit(4, "mm") * nrow(scaled.data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
             show_row_names = TRUE,
             heatmap_legend_param = list(direction = "horizontal")
),
heatmap_legend_side = "bottom")

# dev.off()

