library(tidyverse)
#library("scales")    
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(fs)
library(WriteXLS)
library(compositions)
library(ggpattern)

# load cafog corrected data -----------------------------------------------

abundance_data <- NULL
# Specify the directory path
directory_path <- "analysis/cafog/"

# List all directories in the specified path
folders <- dir_ls(directory_path, type = "directory")

for (folder in folders) {
  file_path <- paste0(folder,"/results.csv")
  
  abundance_data <- rbind(abundance_data,
                          read_csv(file_path,
                                   n_max = 10) %>%
                            mutate(experiment_tp = str_extract(folder, "[^/]+$")) %>%
                            {.}
  )
}

corr_abundance_data <- abundance_data %>%
  separate(
    glycoform,
    into = c("glycoform1", "glycoform2",  "glycoform3"),
    sep = "\\s+or\\s+",
    remove = FALSE
  ) %>%
  select(glycoform1, corr_abundance, corr_abundance_error, experiment_tp) %>%
  mutate(glycoform1 = str_replace_all(glycoform1, "A2", ""))  %>%
  mutate(glycoform1 = str_replace_all(glycoform1, c("G0F/G2F" = "G1F/G1F", 
                                                    "G2F/none" = "none/G2F",
                                                    "G1F/none" = "none/G1F", 
                                                    "G0F/none" = "none/G0F", 
                                                    "G0/none" = "none/G0")
                                      )) %>%
  mutate(glycoform1 = factor(glycoform1, levels = c("none/G0F",
                                                    "none/G1F",
                                                    "none/G2F",
                                                    "G0/G0",
                                                    "G0/G0F",
                                                    "G0F/G0F",
                                                    "G0F/G1F",
                                                    "G1F/G1F",
                                                    "G1F/G2F",
                                                    "G2F/G2F")
                             )) %>%
  separate(experiment_tp,
           into = c("experiment",
                    "timepoint"
           ),
           sep = "_",
           remove = FALSE) %>%
  mutate(experiment = factor(experiment, levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20"))) %>%
  # filter(experiment != "E17")
  # drop_na()
  {.}

save(corr_abundance_data, file = "analysis/corr_abundance_data.RData")

load("analysis/corr_abundance_data.RData")

correct_order <- c("none/G0F",
                   "none/G1F",
                   "none/G2F",
                   "G0/G0",
                   "G0/G0F",
                   "G0F/G0F",
                   "G0F/G1F",
                   "G1F/G1F",
                   "G1F/G2F",
                   "G2F/G2F")



corr_abundance_data <- corr_abundance_data %>%
  # drop_na() %>%
  mutate(experiment = factor(experiment, 
                             levels = c("E13", "E15", "E17", "E19", "E14", "E16", "E18", "E20")
                             )) %>%
  # filter(experiment != c("E17")) %>%
  group_by(experiment, timepoint) %>%
  arrange(match(glycoform1, correct_order)) %>%
  ungroup() %>%
  arrange(experiment_tp) %>%
  {.}

write_csv(x = corr_abundance_data,
          file = "analysis/corr_abundance_data_4tp.csv")
WriteXLS(x = corr_abundance_data,
         ExcelFileName = "analysis/corr_abundance_data_4tp.xls")


# prepare data for differential analysis --------------------------------------------------
# data.matrix <- readxl::read_excel("data/data_matrix_TR.xlsx") %>%
#   column_to_rownames('...1') %>%
#   as.matrix() 
  
data.matrix <- corr_abundance_data %>%
  select(glycoform1, corr_abundance, experiment_tp) %>%
  pivot_wider(values_from = corr_abundance,
              names_from = experiment_tp) %>%
  column_to_rownames('glycoform1') %>%
  as.matrix() 


# filter 144 and 288 tp (just as a comparison with peptide mapping)
data.matrix <- corr_abundance_data %>%
  filter(timepoint %in% c("144", "288")) %>%
  select(glycoform1, corr_abundance, experiment_tp) %>%
  pivot_wider(values_from = corr_abundance,
              names_from = experiment_tp) %>%
  column_to_rownames('glycoform1') %>%
  as.matrix() 

# Set all negative values to 0
data.matrix[data.matrix < 0] <- 0

# Apply log2 transformation (adding 1 to avoid log2(0))
log2_data.matrix <- log2(t(data.matrix + 1)) # should not be transposed, but actually it does not matter! Applies log2 to all numbers and does not care whether samples in rows or in columns. 

#Perform log2 transformation
log2_data.matrix <- t(as.matrix(log2_data.matrix))

log2_data.matrix

meta <- tibble(sample_name = colnames(data.matrix)) %>%
  separate(col = sample_name,
           into = c('experiment', 'timepoint'),
           sep = "_",
           remove = FALSE
  ) %>%
  mutate(condition = case_when(
    experiment %in% c('E13', 'E15', 'E17', 'E19') ~ 'constant',
    experiment %in% c('E14', 'E16', 'E18', 'E20') ~ 'tshifted',
    TRUE ~ 'other'  # This handles any other experiments, if applicable
  ),bioprocess_batch = case_when(
    experiment %in% c('E13', 'E14', 'E15', 'E16') ~ '1',
    experiment %in% c('E17', 'E18', 'E19', 'E20') ~ '2',
    TRUE ~ 'other'  # This handles any other experiments, if applicable
  ))

save(log2_data.matrix, data.matrix, meta, file = "analysis/e13_e20_nglycans_TR.RData")
# log2_data.matrix_TR <- log2_data.matrix
# data.matrix_TR <- data.matrix

# clr transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix

save(clr_data.matrix, clr_data.matrix, meta, file = "analysis/e13_e20_nglycans_clr.RData")


##
data.matrix_tosave <- data.matrix %>% 
  as.data.frame() %>%
  mutate(modcom = rownames(data.matrix)) 

write_csv(data.matrix_tosave, 
          file = "analysis/corr_abundance_data_matrix.csv")
WriteXLS(x = data.matrix_tosave,
         ExcelFileName = "analysis/corr_abundance_data_matrix.xls")

# plot char runs data -----------------------------------------------------
# Define the colors
# color_mapping_experiment <- c(
#   "E13" = "#e41a1c",
#   "E14" = "#377eb8",
#   "E15" = "#4daf4a",
#   "E16" = "#984ea3",
#   "E17" = "#ff7f00",
#   "E18" = "#ffff33",
#   "E19" = "#a65628",
#   "E20" = "#f781bf"
# )

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
# #from Larissa
# c("#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603")
# c("#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F")


plot_bars <- function(data,
                      title = "Fractional abundance",
                      row_number = 1){
  ggplot(data, aes(x = glycoform1, y = corr_abundance, fill = experiment)) +
    geom_col(
      position = position_dodge(width = 0.9)  
    ) + 
    geom_errorbar(
      aes(
        ymin = corr_abundance - corr_abundance_error,
        ymax = corr_abundance + corr_abundance_error,
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
                       limits = c(0,65)) +
    xlab("") +
    theme_bw() +
    theme(text = element_text(size = 16, 
                              face = "bold",
                              family = "sans"),
          axis.text.x = element_text(angle = 90, 
                                     vjust = .5, 
                                     hjust = 1),
          axis.text = element_text(colour = "black"),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          panel.border = element_blank()
          ) +
    ggtitle(title) +
    NULL
  }
  

## plot all experiments
plot_bars(corr_abundance_data,
          title = "Fractional abundance of glycans in all experiments",
          row_number = 2)

ggsave(filename = "figures/corrected_frac_ab_barplot_all_experiments_bold16.png",
       height = 200,
       width = 250,
       units = "mm",
       dpi = 600)

## plot constant t experiments
constant_data <- corr_abundance_data %>%
  filter(experiment %in% c("E13", "E15", "E19"))
plot_bars(constant_data,
          title = "Fractional abundance of glycans in constant temperature experiments",
          row_number = 2)

ggsave(filename = "figures/corrected_frac_ab_barplot_minus17_constant_experiments.png",
       height = 200,
       width = 250,
       units = "mm",
       dpi = 600)

## plot shifted t experiments
tshifted_data <- corr_abundance_data %>%
  filter(experiment %in% c("E14", "E16", "E18", "E20"))
plot_bars(tshifted_data,
          title = "Fractional abundance of glycans in temperature shifted experiments",
          row_number = 2)

ggsave(filename = "figures/corrected_frac_ab_barplot_tshifted_experiments.png",
       height = 200,
       width = 250,
       units = "mm",
       dpi = 600)

## plot all experiments, only 4 timepoints
four_tp_data <-  corr_abundance_data %>% filter(timepoint %in% c("120", "216", "288","336"))

plot_bars(four_tp_data,
          title = "Fractional abundance of glycans in all experiments",
          row_number = 1)

ggsave(filename = "figures/corrected_frac_ab_barplot_minus17_all_experiments_4tp.png",
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)

## plot constant t experiments, only 4 timepoints
constant_four_tp_data <- four_tp_data %>%
  filter(experiment %in% c("E13", "E15", "E19"))
plot_bars(constant_four_tp_data,
          title = "Fractional abundance of glycans in constant temperature experiments",
          row_number = 1)

ggsave(filename = "figures/corrected_frac_ab_barplot_minus17_constant_experiments_4tp.png",
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)

## plot shifted t experiments, only 4 timepoints
tshifted_four_tp_data <- four_tp_data %>%
  filter(experiment %in% c("E14", "E16", "E18", "E20"))
plot_bars(tshifted_four_tp_data,
          title = "Fractional abundance of glycans in temperature shifted experiments",
          row_number = 1)

ggsave(filename = "figures/corrected_frac_ab_barplot_tshifted_experiments_4tp.png",
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)

# plot line plots for each glycan separately -------------------------------

str(corr_abundance_data) # Check the structure of your dataframe
corr_abundance_data$timepoint <- as.numeric(as.character(corr_abundance_data$timepoint))

plot_over_time <- function(data,
                           which_experiment = c("E14")) {
ggplot(data %>% filter(experiment %in% which_experiment), 
       aes(x = timepoint, y = corr_abundance, color = experiment)) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) + # Remove fullrange = TRUE
  scale_color_manual(values = color_mapping_experiment) +
  facet_wrap(~glycoform1, 
             scales = "free_y",
             nrow = 2)
  
  ggsave(filename = paste0("figures/corrected_frac_ab_lineplot_over_time",which_experiment,".png"),
         height = 100,
         width = 250,
         units = "mm",
         dpi = 600)

ggplot(data %>% filter(experiment %in% which_experiment), 
       aes(x = timepoint, y = corr_abundance, fill = experiment)) +
  geom_col(
    position = position_dodge(width = 0.9)  
  ) + 
  geom_errorbar(
    aes(
      ymin = corr_abundance - corr_abundance_error,
      ymax = corr_abundance + corr_abundance_error,
      group = experiment
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_manual(values = color_mapping_experiment, 
                    breaks = names(color_mapping_experiment)) +
  facet_wrap(~glycoform1, 
             scales = "free_y",
             nrow = 2)

ggsave(filename =  paste0("figures/corrected_frac_ab_barplot_over_time",which_experiment,".png"),
       height = 100,
       width = 250,
       units = "mm",
       dpi = 600)
}

plot_over_time(corr_abundance_data, which_experiment = c("E20"))


# plot as heatmap ---------------------------------------------------------
# calculate z-score & plot heatmap 

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



# plot bars horizontally ----------------------------------------------------
subset <- corr_abundance_data %>%
  filter(timepoint %in% c("120", "336")) %>%
  mutate(condition = case_when(
    experiment %in% c("E13", "E15", "E17", "E19") ~ "Constant",
    experiment %in% c("E14", "E16", "E18", "E20") ~ "Temp. shifted",
    TRUE ~ NA_character_  # optional, for any experiments not matched
    )
  ) %>%
  mutate(condition_tp = paste(condition, timepoint, sep = "_")) %>%
  filter(experiment %in% c("E13", "E14", "E15", "E18", "E19", "E20"))
  {}

correct_order_coord_flip <- c("G2F/G2F",
                              "G1F/G2F",
                              "G1F/G1F",
                              "G0F/G1F",
                              "G0F/G0F",
                              "G0/G0F",
                              "G0/G0",
                              "none/G2F",
                              "none/G1F",
                              "none/G0F"
)

color_mapping_condition <- c(
"Constant" = "#E6641E",
"Temp. shifted" = "#4B288C"
)

color_mapping_condition_tp <- c(
  "Constant_120" = "#E6641E",
  "Constant_336" = "#E6641E",
  "Temp. shifted_120" = "#4B288C",
  "Temp. shifted_336" = "#4B288C"
  
)

subset_stats <- subset %>%
  group_by(glycoform1, condition, timepoint) %>%
  summarise(
    mean_frac_ab = mean(corr_abundance),
    sd_frac_ab = sd(corr_abundance)
  ) %>%
  ungroup() %>%
  mutate(
    glycoform1 = factor(glycoform1, 
                        levels = correct_order_coord_flip),
    pattern = ifelse(str_detect(timepoint, "336"), "none", "stripe"),
    condition_tp = paste(condition, timepoint, sep = "_")
  ) %>%
  mutate(pattern = as.character(pattern))

dodge <- position_dodge(width = 0.9)


ggplot(subset_stats) +
  geom_col_pattern(
    aes(
      x = glycoform1,
      y = mean_frac_ab, 
      fill = condition,               
      group = condition_tp,              
      pattern = timepoint,
      pattern_density = timepoint
      ),
    color = "black",
    position = dodge,
    linewidth = 0.025,
    pattern_fill = alpha("black", 0.5),  # Add transparency to the stripes
    # pattern_fill = "black",  
    pattern = "stripe",
    pattern_spacing = 0.025,
    pattern_angle = 45,
    # pattern_density = 0.1
  ) +
  # geom_point(
  #   data = subset,
  #   aes(x = glycoform1,
  #       y = corr_abundance,
  #       group = condition_tp,
  #       fill = condition,
  #       ),
  #   position = dodge,
  #   color = "black",
  #   size = 2,
  #   shape = 21
  # ) +
  geom_errorbar(
    data = subset_stats,
    aes(
      x = glycoform1,
      ymin = mean_frac_ab - sd_frac_ab,
      ymax = mean_frac_ab + sd_frac_ab,
      group = condition_tp
    ),
    position = dodge,
    width = 0.5,
    linewidth = 0.25
  ) +
  scale_fill_manual(values = color_mapping_condition) +
  xlab("") +
  ylab("fractional abundance (%)") +
  labs(title = "Hexose bias corrected glycoforms - intact") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 12,
                            face = "bold",
                            family = "sans"),
        axis.text.y = element_text(colour = "black", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
        ) +
  guides(
    fill = guide_legend(title = "Condition",
                        override.aes = list(pattern = "none")),  
    pattern = guide_legend(title = "Timepoint",
                          override.aes = list(fill = "white", color = "black"))  
)
  
  
ggsave("figures/corrected_frac_ab_tp_120_336_vertical_mean_sd_minusE17E16.png",
  width = 170,
  height = 150,
  units = "mm",
  dpi = 600,
  bg = "transparent"
  )


# no coord flip -----------------------------------------------------------
subset <- corr_abundance_data %>%
  filter(timepoint %in% c("120", "336")) %>%
  mutate(condition = case_when(
    experiment %in% c("E13", "E15", "E17", "E19") ~ "Constant",
    experiment %in% c("E14", "E16", "E18", "E20") ~ "Temp. shifted",
    TRUE ~ NA_character_  # optional, for any experiments not matched
  )
  ) %>%
  mutate(condition_tp = paste(condition, timepoint, sep = "_")) %>%
  filter(experiment %in% c("E13", "E14", "E15", "E18", "E19", "E20"))

correct_order <- c("none/G0F",
                   "none/G1F",
                   "none/G2F",
                   "G0/G0",
                   "G0/G0F",
                   "G0F/G0F",
                   "G0F/G1F",
                   "G1F/G1F",
                   "G1F/G2F",
                   "G2F/G2F")


color_mapping_condition <- c(
  "Constant" = "#E6641E",
  "Temp. shifted" = "#4B288C"
)

color_mapping_condition_tp <- c(
  "Constant_120" = "#E6641E",
  "Constant_336" = "#E6641E",
  "Temp. shifted_120" = "#4B288C",
  "Temp. shifted_336" = "#4B288C"
  
)

subset_stats <- subset %>%
  group_by(glycoform1, condition, timepoint) %>%
  summarise(
    mean_frac_ab = mean(corr_abundance),
    sd_frac_ab = sd(corr_abundance)
  ) %>%
  ungroup() %>%
  mutate(
    glycoform1 = factor(glycoform1, 
                        levels = correct_order),
    pattern = ifelse(str_detect(timepoint, "336"), "none", "stripe"),
    condition_tp = paste(condition, timepoint, sep = "_")
  ) %>%
  mutate(pattern = as.character(pattern))

dodge <- position_dodge(width = 0.9)


ggplot(subset_stats) +
  geom_col_pattern(
    aes(
      x = glycoform1,
      y = mean_frac_ab, 
      fill = condition,               
      group = condition_tp,              
      pattern = timepoint,
      pattern_density = timepoint
    ),
    # color = "black",
    position = dodge,
    linewidth = 0.025,
    pattern_fill = alpha("black", 0.5),  # Add transparency to the stripes
    # pattern_fill = "black",  
    pattern = "stripe",
    pattern_spacing = 0.025,
    pattern_angle = 45,
    # pattern_density = 0.1
  ) +
  # geom_point(
  #   data = subset,
  #   aes(x = glycoform1,
  #       y = corr_abundance,
  #       group = condition_tp,
  #       fill = condition,
  #   ),
  #   position = dodge,
  #   color = "black",
  #   size = 2,
  #   shape = 21
  # ) +
  geom_errorbar(
    data = subset_stats,
    aes(
      x = glycoform1,
      ymin = mean_frac_ab - sd_frac_ab,
      ymax = mean_frac_ab + sd_frac_ab,
      group = condition_tp
    ),
    position = dodge,
    width = 0.5,
    linewidth = 0.25
  ) +
  scale_fill_manual(values = color_mapping_condition) +
  xlab("") +
  ylab("fractional abundance (%)") +
  labs(title = "Hexose bias corrected glycoforms - intact") +
  geom_hline(yintercept = 0, linewidth = .35) +
  # coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 12,
                            face = "bold",
                            family = "sans"),
        axis.text.y = element_text(colour = "black", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   vjust = .5, 
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  ) +
  guides(
    fill = guide_legend(title = "Condition",
                        override.aes = list(pattern = "none")),  
    pattern = guide_legend(title = "Timepoint",
                           override.aes = list(fill = "white", color = "black"))  
  )


ggsave("figures/corrected_frac_ab_tp_120_336_horizontal_mean_sd_minusE17E16.png",
       width = 150,
       height = 150,
       units = "mm",
       dpi = 600,
       bg = "transparent"
)

