library(SplineOmics)
library(readxl) # for loading Excel files
library(here) # For managing filepaths
library(dplyr) # For data manipulation


# loading files -----------------------------------------------------------


load("analysis/e13_e20_nglycans_clr.RData")

# data <- readRDS(xzfile(system.file(
#   "extdata",
#   "proteomics_data.rds.xz",
#   package = "SplineOmics"
# )))
data <- clr_data.matrix


# meta <- read_excel(
#   system.file(
#     "extdata",
#     "proteomics_meta.xlsx",
#     package = "SplineOmics"
#   )
# )

# # Extract the annotation part from the dataframe.
# first_na_col <- which(is.na(data[1, ]))[1]
# annotation <- data |>
#   dplyr::select((first_na_col + 1):ncol(data)) |>
#   dplyr::slice(-c(1:3))
annotation <- rownames(data)

print(head(data))
print(meta)
print(annotation)



# report info -------------------------------------------------------------

report_info <- list(
  omics_data_type = "clr transformed fractinal abundance nglycan data",
  data_description = "nglycan data in BR_24_E13-20",
  data_collection_date = "September 2024",
  analyst_name = "Veronika Schäpertöns",
  contact_info = "veronika.schaepertoensr@plus.ac.at",
  project_name = "DGTX"
)

report_dir <- here::here(
  "analysis",
  "splineomics",
  "explore_data"
)


# Splineomics object ------------------------------------------------------

# splineomics now contains the SplineOmics object.
splineomics <- SplineOmics::create_splineomics(
  data = data,
  meta = meta,
  annotation = annotation,
  report_info = report_info,
  condition = "condition", # Column of meta that contains the levels.
  meta_batch_column = "experiment", # For batch effect removal
  meta_batch2_column = "bioprocess_batch"
)

# Special print.SplineOmics function leads to selective printing
print(splineomics)


# EDA ---------------------------------------------------------------------

plots <- SplineOmics::explore_data(
  splineomics = splineomics, # SplineOmics object
  report_dir = report_dir
)
# hyperparameters testing ------------------------------------------------


# update Splineomics object -----------------------------------------------
splineomics <- SplineOmics::update_splineomics(
  splineomics = splineomics,
  design = "~ 1 + condition*timepoint + batch", # best design formula
  mode = "integrated", # means limma uses the full data for each condition.
  data = data, # data without "outliers" was better
  meta = meta,
  spline_params = list(
    spline_type = c("n"), # natural cubic splines (take these if unsure)
    dof = c(2L) # If you are unsure about which dof, start with 2 and increase
  )
)


# run splineomics ---------------------------------------------------------

splineomics <- SplineOmics::run_limma_splines(
  splineomics = splineomics
) 
