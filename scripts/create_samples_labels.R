

# Define the ranges for each column
experiments <- paste0("E", 13:20)
time_point_label <- paste0("TP", 1:59)
technical_replicates <- 1:3

# Create the data frame using expand.grid
df <- expand.grid(
  experiment = experiments,
  TP = time_point_label,
  technical_replicate = technical_replicates
)

# Define the new time points
time_points <- c(
  1, 8, 24, 32, 48, 56, 72, 74, 76, 80, 96, 100, 104, 120, 122, 124, 
  128, 144, 148, 152, 168, 170, 172, 176, 192, 196, 200, 216, 
  218, 220, 224, 240, 244, 248, 264, 266, 268, 272, 288, 292, 296, 
  312, 314, 316, 320, 336, 340, 344, 360, 362, 364, 368, 384, 392, 
  408, 416, 432, 440, 456
)

# Create a mapping between time points and new time points
time_point_mapping <- data.frame(
  TP = paste0("TP", 1:59),
  time_point = time_points
)

# Merge the mapping with the original data frame
df <- merge(df, time_point_mapping, by = "TP") %>%
  mutate(ymd = '20240606',
         initials = 'TB',
         product = 'cnistcho',
         pngase = 'none') %>%
  select(ymd, initials, product, experiment,TP, time_point, pngase,technical_replicate) %>%
  mutate(sample_name = paste(ymd, initials, product, experiment, time_point, pngase, technical_replicate, sep = "_"))

write_csv(df, "samples_overview.csv")

df_filtered <- df %>%
  filter(time_point %in% c("120", "144", "168", "192", "216", "288", "336")) %>%
  arrange(sample_name, as.numeric(time_point), technical_replicate)

write_csv(df_filtered, "samples_overview_filtered.csv")         
