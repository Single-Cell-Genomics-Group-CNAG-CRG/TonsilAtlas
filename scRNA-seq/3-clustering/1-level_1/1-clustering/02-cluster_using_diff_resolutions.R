# This script clusters the cells in the integrated Seurat object using different
# resolutions.


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_utils <- here::here("scRNA-seq/bin/utils.R")
path_to_tmp <- here::here("scRNA-seq/3-clustering/1-level_1/tmp")
path_to_obj <- here::here("scRNA-seq/results/R_objects/tonsil_rna_integrated_knn.rds")
path_to_save_clusters_df <- str_c(path_to_tmp, "clusters_level_1_df.rds", sep = "/")
path_to_save_harmony_coords <- str_c(path_to_tmp, "harmony_coords.rds", sep = "/")


# Source functions
source(path_to_utils)


# Load data
tonsil <- readRDS(path_to_obj)


# Find clusters at different resolutions
clusters_df <- cluster_diff_resolutions(
  seurat_obj = tonsil,
  min_resolution = 0.01,
  max_resolution = 0.7,
  step = 0.02
)


# Save
dir.create(path_to_tmp, showWarnings = FALSE)
saveRDS(clusters_df, path_to_save_clusters_df)
harmony_coords <- tonsil@reductions$harmony@cell.embeddings[, 1:30]
saveRDS(harmony_coords, path_to_save_harmony_coords)
