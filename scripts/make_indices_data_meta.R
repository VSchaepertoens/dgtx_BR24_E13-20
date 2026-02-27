
library(tidyverse)


load("analysis/charrun_E13-E20_galactosylationindex_V01_20260114_VS.RData")

# Make data table and correct meta and save --------------

#-------------galactosylation_index-------------
#remove E17 288 and 312 timepoints
gi_summary_filt <- gi_summary %>%
  filter(!experiment_tp %in% c("E17_288", "E17_312"))

gi_stats_filt <- gi_summary_filt %>%
  group_by(condition, tp) %>%
  summarise(
    mean_GI = mean(GI),
    sd_GI = sd(GI),
    se_GI = sd(GI, na.rm = TRUE)/sqrt(n()), 
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate (condition_tp = paste(condition, tp, sep = "_"))

#build data table
data_1 <- gi_stats_filt %>%
  select(condition_tp, mean_GI) %>%
  pivot_wider(
    names_from = condition_tp, 
    values_from = mean_GI) %>%
  mutate(feature_id = "galactosylation_index") %>%
  relocate(feature_id)

meta_1 <- tibble(sample_ID = colnames(data_1[,-1])) %>%
  separate(col = sample_ID,
           into = c('condition', 'timepoint'),
           sep = "_",
           remove = FALSE
  )

##--------fucosylation_index--------------

load("analysis/charrun_E13-E20_fucosylationindex_V01_20260114_VS.RData")

#remove E17 288 and 312 timepoints
gi_summary_filt <- gi_summary %>%
  filter(!experiment_tp %in% c("E17_288", "E17_312"))

gi_stats_filt <- gi_summary_filt %>%
  group_by(condition, tp) %>%
  summarise(
    mean_GI = mean(GI),
    sd_GI = sd(GI),
    se_GI = sd(GI, na.rm = TRUE)/sqrt(n()), 
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate (condition_tp = paste(condition, tp, sep = "_"))

#build data table
data_2 <- gi_stats_filt %>%
  select(condition_tp, mean_GI) %>%
  pivot_wider(
    names_from = condition_tp, 
    values_from = mean_GI) %>%
  mutate(feature_id = "fucosylation_index") %>%
  relocate(feature_id)

meta_2 <- tibble(sample_ID = colnames(data_2[,-1])) %>%
  separate(col = sample_ID,
           into = c('condition', 'timepoint'),
           sep = "_",
           remove = FALSE
  )

# join both data and meta into one
identical(meta_1, meta_2)

data <- rbind(data_1,data_2)
meta <- meta_1

annotation <- tibble(feature_ID = data$feature_id)


save(data, meta, file = "analysis/charrun_E13-E20_cqaindices_V01_20260227_VS.RData")  