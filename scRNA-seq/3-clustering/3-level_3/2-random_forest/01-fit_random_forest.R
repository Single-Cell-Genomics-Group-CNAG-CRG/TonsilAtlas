# This script fits a random forest model using each of the clustering performed
# in the previous notebook.


# Load packages
library(tidyverse)
library(randomForest)
library(doParallel)
library(foreach)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]
res <- args[[2]]


# Define paths
path_to_tmp <- str_c(
  here::here("scRNA-seq/3-clustering/3-level_3/tmp"),
  cell_type,
  sep = "/"
)
path_to_clusters_df <- str_c(
  path_to_tmp,
  "clusters_level_3_df.rds",
  sep = "/"
)
path_to_harmony_coords <- str_c(path_to_tmp, "harmony_coords.rds", sep = "/")
path_to_utils <- here::here("scRNA-seq/bin/utils.R")
path_to_oob_dir <- str_c(path_to_tmp, "oob", sep = "/")
path_to_save_oob_df <- str_c(path_to_oob_dir, "/", "oob_accuracy_res", res, ".rds", sep = "")


# Create directories (if they do not exist)
dir.create(path_to_oob_dir, showWarnings = FALSE)


# Source functions
source(path_to_utils)


# Load data
clusters_df <- readRDS(path_to_clusters_df)
harmony_coords <- readRDS(path_to_harmony_coords)


# Fit random forest model
large_cell_types <- c("NBC_MBC", "GCBC", "CD4_T")
n_cells_per_class <- ifelse(cell_type %in% large_cell_types, 3000, 1000)
print("Number of cells per class:")
print(n_cells_per_class)
oob_accuracy <- run_random_forest(
  resolutions_df = clusters_df,
  resolution = res,
  dim_red_coords = harmony_coords,
  n_cells_per_class = n_cells_per_class,
  seed = 123,
  n_trees = 600,
  n_cpu = 4,
  return  = "oob"
)


# Save
dir.create(str_c(path_to_tmp, "oob", sep = "/"))
oob_df <- data.frame(resolution = res, oob_accuracy = oob_accuracy)
saveRDS(oob_df, path_to_save_oob_df)

