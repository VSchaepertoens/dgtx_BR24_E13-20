library(tidyverse)
library(fs)
library(readxl)
# Load necessary libraries
library(dplyr)
library(ggplot2)


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

# load FBA data -----------------------------------------------------------

fba_data <- read_csv("fba_results/condition_specific/iCHO1766_biomass_producing/reaction_data_results_icho1766_FBA_pFBA.csv")


# load pFBA data ----------------------------------------------------------

# pfba_data <- read_csv("fba_results/condition_specific/iCHO1766_biomass_producing/reaction_data_results_icho1766_pFBA.csv")

# join fva, fba and pFBA into one df


# plot fva, fba and pFBA results of particular reactions --------------

glycolysis_rxn <- c("HEX1", "PFK", "PGI", "GAPD", "TPI")

igg_rxn <- c("igg_hc", "igg_lc", "igg_formation", "DM_igg_g_")


filt_fva_data <- combined_fva_data %>%
  filter(
         Reaction %in% glycolysis_rxn)

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


ggplot(joined_fva_fba_pfba, aes(x = Reaction, y = Flux, color = Experiment)) +
  geom_pointrange(aes(ymin = Min_flux, ymax = Max_flux)) +
  facet_grid(Experiment ~ Window) +
  # facet_wrap(~ Window, ncol = 5) +
  scale_color_manual(values = color_mapping_experiment) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(filt_fba_data, aes(x = Reaction)) +
  geom_point(aes(y = FBA_Flux), color = "blue") +
  geom_point(aes(y = pFBA_Flux), color = "red") +
  facet_grid(Experiment ~ Window) +
  # ylim( -0.1, 0.5) +
  # facet_wrap(~ Window, ncol = 5) +
  # scale_color_manual(values = color_mapping_experiment) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))



# PCA --------------------------------------------------------------------
# Assuming your data is in a dataframe called 'df'
# df should have columns: Experiment, Window, Reaction, Flux

# Step 1: Reshape the data to wide format
# Each reaction will be a column, and each row will be a unique combination of Experiment and Window
df_wide <- fba_data %>%
  select(Experiment, Window, Reaction, FBA_Flux) %>%
  # Aggregate duplicates (if any)
  group_by(Experiment, Window, Reaction) %>%
  summarize(Flux = mean(FBA_Flux, na.rm = TRUE), .groups = 'drop') %>%
  # Ensure all combinations are present
  complete(Experiment, Window, Reaction, fill = list(Flux = 0)) %>%
  # Pivot to wide format
  pivot_wider(names_from = Reaction, values_from = Flux) %>%
  # Create unique row identifiers
  mutate(experiment_window = paste(Experiment, Window, sep = "_")) %>%
  # Move to row names
  column_to_rownames("experiment_window")

# Step 2: Remove columns (reactions) with only zeros
df_wide <- df_wide %>%
  select(where(~ any(. != 0)))  # Keep only columns with at least one non-zero value

# Perform PCA
pca_result <- prcomp(x = df_wide %>% select(-Experiment, -Window),
                     scale. = TRUE)

#  Extract PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$experiment_window <- rownames(df_wide)
pca_scores$Experiment <- df_wide$Experiment
pca_scores$Window <- df_wide$Window

# Step 1: Calculate the variance explained by each PC
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Step 2: Create a summary table
variance_summary <- data.frame(
  PC = paste0("PC", 1:length(variance_explained)),
  Variance_Explained = variance_explained
)

# Display the variance explained by PC1 and PC2
cat("Variance explained by PC1:", round(variance_summary$Variance_Explained[1], 2), "%\n")
cat("Variance explained by PC2:", round(variance_summary$Variance_Explained[2], 2), "%\n")
cat("Variance explained by PC3:", round(variance_summary$Variance_Explained[3], 2), "%\n")
cat("Variance explained by PC4:", round(variance_summary$Variance_Explained[4], 2), "%\n")

# Step 4: Visualize the PCA results
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
pca_per_window(pca_scores, 2)
