library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# loading data ------------------------------------------------------------

constant <-  read_csv("fba_results/nottshifted/case_2.1/reaction_data_results.csv") %>%
  mutate(Condition = "constant")

tshifted <-  read_csv("fba_results/tshifted/case_2.1/reaction_data_results.csv") %>%
  mutate(Condition = "tshifted")

merged_data <- constant %>%
  rbind(tshifted)

# evaluate fluxes by PCA --------------------------------------------------
#first need to create averages
averaged_merged_data <- merged_data %>%
  group_by(Experiment, Window, Reaction, Condition) %>%
  summarise(mean_flux = mean(Flux),
            sd_flux = sd(Flux)) %>%
  ungroup() %>%
  mutate(Experiment_window = paste(Experiment, Window, sep = "_"))
   

data.matrix <- averaged_merged_data %>%
  select(Experiment_window, mean_flux, Reaction) %>%
  pivot_wider(values_from = mean_flux,
              names_from = Experiment_window) %>%
  column_to_rownames('Reaction') %>%
  filter(rowSums(. == 0) != ncol(.)) %>%
  as.matrix()

boxplot(data.matrix,
        las = 2)
#checking the size of the matrix
nrow(data.matrix) 
ncol(data.matrix)

# 1. calculate the principal components
#calculate principal components
results <- prcomp(data.matrix[,1:2], scale = TRUE)

#reverse the signs
results$rotation <- -1*results$rotation

#display principal components
results$rotation


#reverse the signs of the scores
results$x <- -1*results$x

#display the first six scores
head(results$x)

biplot(results, scale = 0)

#calculate total variance explained by each principal component
results$sdev^2 / sum(results$sdev^2)

#calculate total variance explained by each principal component
var_explained = results$sdev^2 / sum(results$sdev^2)

#create scree plot
qplot(c(1:2), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
# subset key core reactions -----------------------------------------------

subset <- merged_data %>%
  filter(Reaction %in% c("HEX1", "PFK", "PGI", "GAPD", "TPI")) %>%
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
f1 = colorRamp2(seq(0,max(abs(data.matrix_subset)),length = 9),
                c("black",
                  "#ffeda0",
                  "#fed976",
                  "#feb24c",
                  "#fd8d3c",
                  "#fc4e2a",
                  "#e31a1c",
                  "#bd0026",
                  "#800026"),
                space = "RGB")

ht = Heatmap(data.matrix_subset,
             col = f1,
             cluster_columns = FALSE,
             rect_gp = gpar(col = "white", lwd = 2),
             name = "flux (mMol/gDCW/hr-1)",
             row_gap = unit(4, "pt"),
             column_gap = unit(4, "pt"),
             # width = unit(4, "mm") * ncol(data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
             # height = unit(4, "mm") * nrow(data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
             show_row_names = TRUE,
             heatmap_legend_param = list(direction = "horizontal"))

draw(ht)





