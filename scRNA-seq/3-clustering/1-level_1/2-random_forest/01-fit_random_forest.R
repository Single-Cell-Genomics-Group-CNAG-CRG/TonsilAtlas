# This script fits a random forest model using each of the clustering performed
# in the previous notebook.


# Load packages
library(tidyverse)
library(randomForest)
library(doParallel)
library(foreach)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
res <- args[[1]]


# Define paths
path_to_utils <- here::here("scRNA-seq/bin/utils.R")
path_to_tmp <- here::here("scRNA-seq/3-clustering/1-level_1/tmp")
path_to_harmony_coords <- str_c(path_to_tmp, "harmony_coords.rds", sep = "/")
path_to_clusters_df <- str_c(path_to_tmp, "clusters_level_1_df.rds", sep = "/")
path_to_save_oob_df <- str_c(path_to_tmp, "/oob/", "oob_accuracy_res", res, ".rds", sep = "")


# Source functions
source(path_to_utils)


# Load data
clusters_df <- readRDS(path_to_clusters_df)
harmony_coords <- readRDS(path_to_harmony_coords)


# Fit random forest model
oob_accuracy <- run_random_forest(
  resolutions_df = clusters_df,
  resolution = res,
  dim_red_coords = harmony_coords,
  n_cells_per_class = 3000,
  seed = 123,
  n_trees = 600,
  n_cpu = 4,
  return  = "oob"
)


# Save
dir.create(str_c(path_to_tmp, "oob", sep = "/"))
oob_df <- data.frame(resolution = res, oob_accuracy = oob_accuracy)
saveRDS(oob_df, path_to_save_oob_df)

