# This script clusters the cells in the integrated Seurat object using different
# resolutions.


# Load packages
library(Seurat)
library(tidyverse)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]


# Define paths
path_to_dir <- str_c(
  here::here("scRNA-seq/results/R_objects/level_3/"),
  cell_type,
  sep = ""
)
path_to_obj <- str_c(
  path_to_dir,
  "/",
  cell_type,
  "_integrated_level_3.rds",
  sep = ""
)
path_to_utils <- here::here("scRNA-seq/bin/utils.R")
path_to_tmp <- here::here("scRNA-seq/3-clustering/3-level_3/tmp")
path_to_tmp_cell_type <- str_c(path_to_tmp, cell_type, sep = "/")
path_to_save_clusters_df <- str_c(
  path_to_tmp_cell_type,
  "clusters_level_3_df.rds",
  sep = "/"
)
path_to_save_harmony_coords <- str_c(
  path_to_tmp_cell_type,
  "harmony_coords.rds",
  sep = "/"
)


# Create directories (if they do not exist)
dir.create(path_to_tmp, showWarnings = FALSE)
dir.create(path_to_tmp_cell_type, showWarnings = FALSE)


# Source functions
source(path_to_utils)


# Load data
seurat <- readRDS(path_to_obj)


# Find clusters at different resolutions
clusters_df <- cluster_diff_resolutions(
  seurat_obj = seurat,
  min_resolution = 0.01,
  max_resolution = 0.7,
  step = 0.02
)


# Save
dir.create(path_to_tmp, showWarnings = FALSE)
saveRDS(clusters_df, path_to_save_clusters_df)
harmony_coords <- seurat@reductions$harmony@cell.embeddings[, 1:30]
saveRDS(harmony_coords, path_to_save_harmony_coords)
