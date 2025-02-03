library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(limma)

# loading data ------------------------------------------------------------
#iCHO1766
constant <-  read_csv("fba_results/nottshifted/case_2.1/reaction_data_results.csv") %>%
  mutate(Condition = "constant",
         gem = "icho1766")

tshifted <-  read_csv("fba_results/tshifted/case_2.1/reaction_data_results.csv") %>%
  mutate(Condition = "tshifted",
         gem = "icho1766")

merged_data <- constant %>%
  rbind(tshifted)

#iCHO2441
constant_icho2441 <-  read_csv("fba_results/nottshifted/case_2.1/reaction_data_results_icho2441.csv") %>%
  mutate(Condition = "constant",
         gem = "icho2441")

tshifted_icho2441 <-  read_csv("fba_results/tshifted/case_2.1/reaction_data_results_icho2441.csv") %>%
  mutate(Condition = "tshifted",
         gem = "icho2441")

merged_data_icho2441 <- constant_icho2441 %>%
  rbind(tshifted_icho2441)
# evaluate fluxes by PCA --------------------------------------------------
#first need to create averages
averaged_merged_data <- merged_data %>%
  group_by(Experiment, Window, Reaction, Condition, gem) %>%
  summarise(mean_flux = mean(Flux),
            sd_flux = sd(Flux)) %>%
  ungroup() %>%
  mutate(Experiment_window = paste(Experiment, Window, gem, sep = "_"))

averaged_merged_data_icho2441 <- merged_data_icho2441 %>%
  group_by(Experiment, Window, Reaction, Condition, gem) %>%
  summarise(mean_flux = mean(Flux),
            sd_flux = sd(Flux)) %>%
  ungroup() %>%
  mutate(Experiment_window = paste(Experiment, Window, gem,  sep = "_"))
   

data.matrix <- averaged_merged_data %>%
  select(Experiment_window, mean_flux, Reaction) %>%
  pivot_wider(values_from = mean_flux,
              names_from = Experiment_window) %>%
  column_to_rownames('Reaction') %>%
  # filter(rowSums(. == 0) != ncol(.)) %>%
  as.matrix()

data.matrix_icho2441 <- averaged_merged_data_icho2441 %>%
  select(Experiment_window, mean_flux, Reaction) %>%
  pivot_wider(values_from = mean_flux,
              names_from = Experiment_window) %>%
  column_to_rownames('Reaction') %>%
  # filter(rowSums(. == 0) != ncol(.)) %>%
  as.matrix()

max(data.matrix)
min(data.matrix)

max(data.matrix_icho2441)
min(data.matrix_icho2441)

graphics::boxplot(data.matrix, las = 2)
graphics::boxplot(data.matrix_icho2441, las = 2)
#checking the size of the matrix
nrow(data.matrix) 
ncol(data.matrix)
nrow(data.matrix_icho2441) 
ncol(data.matrix_icho2441)

plotMDS(data.matrix,
        gene.selection = "common",
        var.explained = "TRUE")

plotMDS(data.matrix_icho2441,
        gene.selection = "common",
        var.explained = "TRUE")


# matching dataframes by Reaction and merging -----------------------------
data_df1 <- as.data.frame(data.matrix) %>%
  rownames_to_column(var = "Reaction") %>%
  # filter(rowSums(. == 0) != ncol(.)) %>%
  {.}
  
data_df2 <- as.data.frame(data.matrix_icho2441) %>%
  rownames_to_column(var = "Reaction") %>%
  # filter(rowSums(. == 0) != ncol(.)) %>%
  {.}
  

merged_data_df <- merge(data_df1, data_df2, by = "Reaction")

merged_matrix <- merged_data_df %>%
  column_to_rownames("Reaction") %>%
  filter(rowSums(. == 0) != ncol(.)) %>%
  as.matrix()
nrow(merged_matrix) 

meta <- data.frame(sample_name = colnames(merged_matrix)) %>%
  separate(sample_name, c("Experiment", "Window", "gem"), sep = "_", remove = FALSE)
  
save(merged_matrix, meta, file = "fba_results/merged_fluxes.R")

nrow(merged_matrix) 
ncol(merged_matrix)

meta_w1 <- meta %>%
filter(Window == 1)

meta_w2 <- meta %>%
  filter(Window == 2)

meta_w3 <- meta %>%
  filter(Window == 3)

meta_w4 <- meta %>%
  filter(Window == 4)

meta_w5 <- meta %>%
  filter(Window == 5)


plotMDS(merged_matrix[,meta_w1$sample_name],
        gene.selection = "common",
        var.explained = "TRUE")

plotMDS(merged_matrix[,meta_w2$sample_name],
        gene.selection = "common",
        var.explained = "TRUE")

plotMDS(merged_matrix[,meta_w3$sample_name],
        gene.selection = "common",
        var.explained = "TRUE")

plotMDS(merged_matrix[,meta_w4$sample_name],
        gene.selection = "common",
        var.explained = "TRUE")

plotMDS(merged_matrix[,meta_w5$sample_name],
        gene.selection = "common",
        var.explained = "TRUE")


# plot pearson correlations -----------------------------------------------


# subset key core reactions -----------------------------------------------

subset <- merged_data %>%
  filter(Reaction %in% c("HEX1", "PFK", "PGI", "GAPD", "TPI")) %>%
  group_by(Experiment, Window, Reaction, Condition) %>%
  summarise(mean_flux = mean(Flux),
            sd_flux = sd(Flux)) %>%
  ungroup() %>%
  mutate(Experiment_window = paste(Experiment, Window, sep = "_"))

subset <- merged_data_icho2441 %>%
  # filter(Reaction %in% c("igg_hc", "igg_lc", "igg_formation","DM_igg_g_")) %>%
  filter(Reaction %in% c("igg_hc", "igg_lc", "igg_formation", "DM_igg[g]")) %>%
  group_by(Experiment, Window, Reaction, Condition) %>%
  summarise(mean_flux = mean(Flux),
            sd_flux = sd(Flux)) %>%
  ungroup() %>%
  mutate(Experiment_window = paste(Experiment, Window, sep = "_"))


data.matrix_subset <- subset %>%
  select(Experiment_window, mean_flux, Reaction) %>%
  pivot_wider(values_from = mean_flux,
              names_from = Experiment_window) %>%
  column_to_rownames('Reaction') %>%
  as.matrix()

arrange(c(2,4,3,5,1))



# plot heatmap ------------------------------------------------------------
#heatmap settings
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
#color scheme
f1 = colorRamp2(seq(-max(abs(data.matrix_subset)),max(abs(data.matrix_subset)),length = 9),
                c("#4575b4",
                  "#74add1",
                  "#abd9e9",
                  "#e0f3f8",
                  "grey",
                  "#fee090",
                  "#fdae61",
                  "#f46d43",
                  "#d73027"),
                space = "RGB")

png(filename = "figures/fba/intracellular_reactions_igg_formation_icho2441_biomass.png",
    width = 300,
    height = 50,
    units = "mm",
    res = 300)

ht = Heatmap(data.matrix_subset,
             col = f1,
             cluster_columns = FALSE,
             rect_gp = gpar(col = "white", lwd = 2),
             name = "flux (mMol/gDCW/hr-1)",
             row_gap = unit(4, "pt"),
             column_gap = unit(4, "pt"),
             width = unit(4, "mm") * ncol(data.matrix_subset) + 5 * unit(4, "pt"), # to make each cell a square
             height = unit(4, "mm") * nrow(data.matrix_subset) + 5 * unit(4, "pt"), # to make each cell a square
             show_row_names = TRUE,
             heatmap_legend_param = list(direction = "horizontal"))


draw(ht)
dev.off()


# plot heatmap of the whole dataset ---------------------------------------
nrow(merged_matrix)

# Logical matrix indicating TRUE for values smaller than -0.05 or larger than 0.05
out_of_range_matrix <- data.matrix < -0.05 | data.matrix > 0.05

# Keep rows that contain at least one value satisfying the condition
rows_to_keep <- rowSums(out_of_range_matrix) > 0
filtered_matrix <- data.matrix[rows_to_keep, ]

# Keep columns that contain at least one value satisfying the condition
cols_to_keep <- colSums(out_of_range_matrix) > 0
filtered_matrix <- filtered_matrix[, cols_to_keep]

# Display filtered matrix
filtered_matrix




# Display the resulting matrix
matrix_005

ht = Heatmap(filtered_matrix, 
             show_row_names = FALSE,
             name = "flux (mMol/gDCW/hr-1)",
             cluster_columns = FALSE)
draw(ht)

max(filtered_matrix)
min(filtered_matrix)


