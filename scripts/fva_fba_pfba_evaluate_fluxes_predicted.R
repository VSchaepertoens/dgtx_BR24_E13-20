library(tidyverse)
library(fs)
library(readxl)
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)


# load FVA data and combine into a single df ------------------------------
df <- tibble(fva_full_path = dir_ls(path = "fba_results/condition_specific/fva_iCHO1766/",regexp =  ".*\\.xlsx"),) %>%
  separate(fva_full_path,
           into = c("dir1", "dir2", "dir3", "filename"),
           sep = "/",
           remove = FALSE) %>%
  separate(filename,
           into = c("method", "experiment_name", "Experiment", "window_name", "Window", "iteration_name", "iteration"),
           sep = "_",
           remove = FALSE) %>%
  select(fva_full_path,filename, Experiment, Window)


# Initialize an empty list to store the processed fva_data tables
fva_data_list <- list()

for (i in 1:nrow(df)) {
  # Get the fva path and additional columns for the current row
  fva_path <- df$fva_full_path[i]
  col1_value <- df$Filename[i]
  col2_value <- df$Experiment[i]
  col3_value <- df$Window[i]
  
  # Load the fva_data from the Excel file
  fva_data <- read_excel(fva_path)
  
  # Add the three columns from df to fva_data
  fva_data <- fva_data %>%
    mutate(Filename = col1_value,
           Experiment = col2_value,
           Window = col3_value)
  
  # Append the processed fva_data to the list
  fva_data_list[[i]] <- fva_data
}

# Combine all fva_data tables into one 
combined_fva_data <- bind_rows(fva_data_list) %>%
  mutate(Window = as.numeric(Window)) %>%
  rename(Min_flux = `Min. Flux`,
         Max_flux = `Max. Flux`)

# load FBA and pFBA data -----------------------------------------------------------

fba_data <- read_csv("fba_results/condition_specific/iCHO1766_biomass_producing/reaction_data_results_icho1766_FBA_pFBA.csv")

# plot fva, fba and pFBA results of particular reactions --------------

glycolysis_rxn <- c("HEX1", "PFK", "PGI", "GAPD", "TPI")

igg_rxn <- c("igg_hc", "igg_lc", "igg_formation", "DM_igg_g_")


filt_fva_data <- combined_fva_data %>%
  filter(
         Reaction %in% igg_rxn)

filt_fba_data <- fba_data %>%
  filter(
         Reaction %in% glycolysis_rxn)

# filt_pfba_data <- pfba_data %>%
#   filter(
#     Reaction %in% glycolysis_rxn) %>%
#   select(Experiment, Window, Reaction, Flux) %>%
#   rename("pfba_flux" = "Flux")
# 
# joined_fva_fba <- left_join(x = filt_fba_data, y = filt_fva_data, by = c("Experiment", "Reaction", "Window"))
# joined_fva_fba_pfba <- left_join(x = joined_fva_fba, y = filt_pfba_data, by = c("Experiment", "Reaction", "Window"))


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
# 
# ggplot(joined_fva_fba_pfba, aes(x = Reaction, y = Flux, color = Experiment)) +
#   geom_pointrange(aes(ymin = Min_flux, ymax = Max_flux)) +
#   facet_grid(Experiment ~ Window) +
#   # facet_wrap(~ Window, ncol = 5) +
#   scale_color_manual(values = color_mapping_experiment) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))

ggplot(filt_fba_data, aes(x = Reaction)) +
  geom_point(aes(y = FBA_Flux), color = "blue") +
  geom_point(aes(y = pFBA_Flux), color = "red") +
  facet_grid(Experiment ~ Window) +
  # ylim( -0.1, 0.5) +
  # facet_wrap(~ Window, ncol = 5) +
  # scale_color_manual(values = color_mapping_experiment) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


# reshaping data ----------------------------------------------------------
# Step 1: Reshape the data to wide format
# Each reaction will be a column, and each row will be a unique combination of Experiment and Window
df_wide_pFBA <- fba_data %>%
  select(Experiment, Window, Reaction, pFBA_Flux) %>%
  # Aggregate duplicates (if any)
  group_by(Experiment, Window, Reaction) %>%
  summarize(Flux = mean(pFBA_Flux, na.rm = TRUE), .groups = 'drop') %>%
  # Ensure all combinations are present
  complete(Experiment, Window, Reaction, fill = list(Flux = 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = Reaction, values_from = Flux) %>%
  # Add column to indicate method type
  mutate(Method = "pfba") %>%
  # Create unique row identifiers
  mutate(experiment_window_method = paste(Experiment, Window, Method,  sep = "_")) %>%
  # Move to row names
  column_to_rownames("experiment_window_method")

df_wide_FBA <- fba_data %>%
  select(Experiment, Window, Reaction, FBA_Flux) %>%
  # Aggregate duplicates (if any)
  group_by(Experiment, Window, Reaction) %>%
  summarize(Flux = mean(FBA_Flux, na.rm = TRUE), .groups = 'drop') %>%
  # Ensure all combinations are present
  complete(Experiment, Window, Reaction, fill = list(Flux = 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = Reaction, values_from = Flux) %>%
  # Add column to indicate method type
  mutate(Method = "fba") %>%
  # Create unique row identifiers
  mutate(experiment_window_method = paste(Experiment, Window, Method,  sep = "_")) %>%
  # Move to row names
  column_to_rownames("experiment_window_method")

# bind both FBA and pFBA
df_wide <- rbind(df_wide_pFBA, df_wide_FBA)

# PCA --------------------------------------------------------------------
# Remove columns (reactions) with only zeros
df_wide <- df_wide %>%
  select(where(~ any(. != 0)))  # Keep only columns with at least one non-zero value

# Perform PCA
pca_result <- prcomp(x = df_wide %>% select(-Experiment, -Window, -Method),
                     scale. = TRUE)

#  Extract PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$experiment_window_method <- rownames(df_wide)
pca_scores$Experiment <- df_wide$Experiment
pca_scores$Window <- df_wide$Window
pca_scores$Method <- df_wide$Method

# Calculate the variance explained by each PC
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create a summary table
variance_summary <- data.frame(
  PC = paste0("PC", 1:length(variance_explained)),
  Variance_Explained = variance_explained
)

# Display the variance explained by PC1 and PC2
cat("Variance explained by PC1:", round(variance_summary$Variance_Explained[1], 2), "%\n")
cat("Variance explained by PC2:", round(variance_summary$Variance_Explained[2], 2), "%\n")
cat("Variance explained by PC3:", round(variance_summary$Variance_Explained[3], 2), "%\n")
cat("Variance explained by PC4:", round(variance_summary$Variance_Explained[4], 2), "%\n")

# Visualize the PCA results
# Plot PCA with color by Experiment
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Experiment)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Flux Data by Experiment",
       x = "Principal Component 1",
       y = "Principal Component 2")

# Plot PCA with color by Window
ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(Window))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Flux Data by Time Window",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Window")

# Plot PCA with color by Window
ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(Method))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Flux Data by Method",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Method")

# # Plot PCA with color by experiment_window
# ggplot(pca_scores, aes(x = PC1, y = PC2, color = experiment_window)) +
#   geom_point(size = 3) +
#   theme_minimal() +
#   labs(title = "PCA of Flux Data by Time Window",
#        x = "Principal Component 1",
#        y = "Principal Component 2",
#        color = "Window")


pca_per_window <- function(scores_from_pca,
                           experimental_window = 1){
# Filter PCA scores for Window 1
pca_scores_window <- scores_from_pca %>% filter(Window == experimental_window)

# Plot PCA for Window 1, colored by Experiment
ggplot(pca_scores_window, aes(x = PC1, y = PC2, color = Experiment)) +
  geom_point(size = 3) +
  scale_color_manual(values = color_mapping_experiment) +
  theme_minimal() +
  labs(title = paste("PCA of Flux Data for Window ",experimental_window),
       x = "Principal Component 1",
       y = "Principal Component 2")
}
pca_per_window(pca_scores, 5)


# Flux correlation analysis -----------------------------------------------
# Analyze the correlation between fluxes from FBA, pFBA, and FVA.
# 
# Steps:
#   
#   Compute the correlation matrix between the fluxes from FBA, pFBA, and FVA (e.g., using Pearson or Spearman correlation).
# 
# Visualize the correlation matrix as a heatmap to identify reactions with highly correlated or uncorrelated fluxes.
# 
# Focus on reactions with high variability (e.g., those with large differences between FVA min and max fluxes).
#pFBA and FBA matrix
fba_pfba_matrix <- df_wide %>%
  select(!c(Experiment, Window, Method)) %>%
  t()

cor_fba_pfba_matrix <- cor(fba_pfba_matrix, method = "spearman")
cor_fba_pfba_matrix[cor_fba_pfba_matrix == 1] <- NA
Heatmap(cor_fba_pfba_matrix, na_col = "grey")

#pFBA matrix
pfba_matrix <- df_wide %>%
  filter(Method == "pfba") %>%
  select(!c(Experiment, Window, Method)) %>%
  t()

cor_pfba_matrix <- cor(pfba_matrix, method = "spearman")
cor_pfba_matrix[cor_pfba_matrix == 1] <- NA
Heatmap(cor_pfba_matrix, na_col = "grey")

#FBA matrix
fba_matrix <- df_wide %>%
  filter(Method == "fba") %>%
  select(!c(Experiment, Window, Method)) %>%
  t()

cor_fba_matrix <- cor(fba_matrix, method = "spearman")
cor_fba_matrix[cor_fba_matrix == 1] <- NA
Heatmap(cor_fba_matrix, na_col = "grey")

# Reaction variability analysis -------------------------------------------
# Use FVA results to identify reactions with high variability (i.e., large differences between min and max fluxes).
# 
# Steps:
#   
#   Calculate the range of variability for each reaction:
#   Variability Range=FVA max−FVA min
# Variability Range=FVA max−FVA min
# 
# Rank reactions by their variability range and focus on the most variable reactions.
# 
# Compare the variability of these reactions across different conditions or experiments.

variability_fva_matrix <- combined_fva_data %>%
  mutate(variability = abs(Max_flux) - abs(Min_flux),
         experiment_window = paste(Experiment, Window, sep = "_")) %>%
  select(Reaction,variability,experiment_window) %>%
  # Pivot to wide format
  pivot_wider(names_from = Reaction, values_from = variability) %>%
  select(where(~ !all(. == 1000))) %>%  # Keep columns where not all values are 1000
  select(where(~ any(. != 0))) %>%     # Keep columns with at least one non-zero value
  column_to_rownames("experiment_window") %>%
  t()

Heatmap(matrix = variability_fva_matrix,
        show_row_names = FALSE,
        cluster_columns = FALSE)


# Flux difference analysis ------------------------------------------------
# Compare the fluxes from FBA and pFBA to identify reactions with significant differences.
# 
# Steps:
#   
#   Compute the absolute difference between FBA and pFBA fluxes for each reaction:
#   Difference=∣FBA flux−pFBA flux∣
# Difference=∣FBA flux−pFBA flux∣
# 
# Rank reactions by their flux differences and focus on those with the largest differences.
# 
# Investigate why these reactions have different fluxes (e.g., due to alternative pathways or degeneracy in the solution space).

# Calculate absolute differences between FBA and pFBA fluxes
  pfba_fba_differences <- fba_data %>%
  mutate(difference = abs(FBA_Flux - pFBA_Flux),
         experiment_window = paste(Experiment, Window, sep = "_")) %>%
  select(Reaction,difference,experiment_window) %>%
  # Pivot to wide format
  pivot_wider(names_from = Reaction, values_from = difference) %>%
  select(where(~ any(. != 0)))  %>% # Keep only columns with at least one non-zero value 
  column_to_rownames("experiment_window") %>%
  t()

Heatmap(matrix = pfba_fba_differences,
        show_row_names = FALSE)



