library(tidyverse)
library(dplyr)
library(janitor)
library(ComplexHeatmap)
library(circlize)

# load data E17 -------------------------------------------------------
# subfolders <- c("E17_originator_TB_KB_trypsin", 
#                 "ptm_E13toE18_originator_TB_KB_trypsin.presets")
subfolder <- c("E17_originator_TB_KB_trypsin")

data_path <- paste0("analysis/peptide_mapping/", 
                    subfolder, 
                    "/results_missed_cleavages.csv")

data <-  read_csv(data_path)

data_filt <- data %>% filter(is.na(`In silico`))

result_E17 <- data_filt %>%
  clean_names() %>%
  group_by(ms_alias_name, digest_missed_cleavages) %>% # Group by sample and missed_cleavage
  summarise(total_val = sum(total_xic_auc_averagine), .groups = "drop") %>% # Count occurrences
  pivot_wider(names_from = digest_missed_cleavages, values_from = total_val, values_fill = 0) %>% # Spread missed_cleavage into columns
  mutate(across(2:4, ~ . / rowSums(across(2:3)) * 100)) %>%
  filter(grepl("E17", ms_alias_name)) %>%
  column_to_rownames("ms_alias_name")
  

# load data E13-E18 -------------------------------------------------------

  subfolder <- c("ptm_E13toE18_originator_TB_KB_trypsin.presets")
  data_path <- paste0("analysis/peptide_mapping/", 
                      subfolder, 
                      "/results_missed_cleavages.csv")
  
  data <-  read_csv(data_path)
  

  data_filt <- data %>% filter(is.na(`In silico`))

  result_E13_E18 <- data_filt %>%
    clean_names() %>%
    group_by(ms_alias_name, digest_missed_cleavages) %>% # Group by sample and missed_cleavage
    summarise(total_val = sum(total_xic_auc_averagine), .groups = "drop") %>% # Count occurrences
    pivot_wider(names_from = digest_missed_cleavages, values_from = total_val, values_fill = 0) %>% # Spread missed_cleavage into columns
    mutate(across(2:4, ~ . / rowSums(across(2:3)) * 100)) %>%
    # filter(grepl("E", ms_alias_name)) %>%
    column_to_rownames("ms_alias_name")

  #just originator
  result_merged <- data_filt %>%
    clean_names() %>%
    group_by(ms_alias_name, digest_missed_cleavages) %>% # Group by sample and missed_cleavage
    summarise(total_val = sum(total_xic_auc_averagine), .groups = "drop") %>% # Count occurrences
    pivot_wider(names_from = digest_missed_cleavages, values_from = total_val, values_fill = 0) %>% # Spread missed_cleavage into columns
    mutate(across(2:4, ~ . / rowSums(across(2:3)) * 100)) %>%
    filter(!(grepl("E", ms_alias_name))) %>%
    # column_to_rownames("ms_alias_name") %>%
    # t() 
  
    write_csv(file = paste0("analysis/peptide_mapping/", 
                            subfolder, 
                            "/results_missed_cleavages_originator.csv"))
    
# load data sva -----------------------------------------------------------
  subfolder <- c("sva57_E13toE17_originator_TB_KB_T_C_G_presets")
  data_path <- paste0("analysis/peptide_mapping/", 
                      subfolder, 
                      "/results_missed_cleavages.csv")
  
  data <-  read_csv(data_path)
  
  data_filt <- data %>% filter(is.na(`In silico`))
  
  result_sva <- data_filt %>%
    clean_names() %>%
    group_by(ms_alias_name, digest_missed_cleavages) %>% # Group by sample and missed_cleavage
    summarise(total_val = sum(total_xic_auc_averagine), .groups = "drop") %>% # Count occurrences
    pivot_wider(names_from = digest_missed_cleavages, values_from = total_val, values_fill = 0) %>% # Spread missed_cleavage into columns
    mutate(across(2:8, ~ . / rowSums(across(2:7)) * 100)) %>%
    # filter(grepl("E", ms_alias_name)) %>%
    column_to_rownames("ms_alias_name")

# merge all together ------------------------------------------------------

  
result_merged <- result_E17 %>%
    rbind(result_E13_E18) %>%
    t()
  
result_merged <-  result_sva %>%
    t()


# plot data as heatmap ----------------------------------------------------
f1 = colorRamp2(seq(0,
                    100,
                    length = 9),
                c("#ffffe5",
                  "#f7fcb9",
                  "#d9f0a3",
                  "#addd8e",
                  "#78c679",
                  "#41ab5d",
                  "#238443",
                  "#006837",
                  "#004529"),
                space = "RGB")
  
  f2 = colorRamp2(seq(0,
                      100,
                      length = 10),
                  c("#8dd3c7",
                    "#ffffb3",
                    "#bebada",
                    "#fb8072",
                    "#80b1d3",
                    "#fdb462",
                    "#b3de69",
                    "#fccde5",
                    "#d9d9d9",
                    "#bc80bd"),
                  space = "RGB")
  

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

png(filename = "figures/peptide_mapping/missed_cleavages_qualitative_SVA.png",
    height = 20,
    width = 40,
    units = "cm",
    res = 600)

draw(Heatmap(result_merged,
             col = f2,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             row_gap = unit(4, "pt"),
             column_gap = unit(4, "pt"),
             width = unit(4, "mm") * ncol(result_merged) + 5 * unit(4, "pt"), # to make each cell a square
             height = unit(4, "mm") * nrow(result_merged) + 5 * unit(4, "pt"),
             heatmap_legend_param = list(
               title = "Percent of peptides with missed cleavages",  # Title for the legend
               title_position = "leftcenter-rot",  # Position the title at the top of the legend
               title_gp = gpar(fontsize = 12, rot = 90),  # Rotate the title by 90 degrees
               at = seq(0, 100, by = 10),  # The annotation values on the legend (0 to 100 in steps of 10)
               labels = seq(0, 100, by = 10)  # Labels corresponding to the annotation values
             ))) # to make each cell a square))

dev.off()


