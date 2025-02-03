library(tidyverse)
#library("scales")    
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(fs)
library(WriteXLS)

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

# Set all negative values to 0
data.matrix[data.matrix < 0] <- 0

# Apply log2 transformation (adding 1 to avoid log2(0))
log2_data.matrix <- log2(t(data.matrix + 1))

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
  ))

save(log2_data.matrix, data.matrix, meta, file = "analysis/e13_e20_nglycans_TR.RData")
# log2_data.matrix_TR <- log2_data.matrix
# data.matrix_TR <- data.matrix

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

