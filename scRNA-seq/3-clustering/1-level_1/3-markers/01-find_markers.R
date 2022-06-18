# This script finds the markers for one cluster of interest (COI) in level 1


# Load packages
library(Seurat)
library(tidyverse)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
coi <- args[[1]]


# Define paths
path_to_obj <- here::here("scRNA-seq/results/R_objects/tonsil_rna_integrated_clustered_level_1.rds")
path_to_save <- here::here("scRNA-seq/3-clustering/1-level_1/tmp/markers_level_1/markers_cluster_")
path_to_save <- str_c(path_to_save, coi, "_level_1.rds", sep = "")


# Load data
tonsil <- readRDS(path_to_obj)
print(tonsil)


# Find Markers
markers_df <- FindMarkers(
  tonsil,
  ident.1 = coi,
  logfc.threshold = 0.4,
  test.use = "wilcox",
  only.pos = TRUE,
  verbose = TRUE
)


# Save Markers
saveRDS(markers_df, path_to_save)
