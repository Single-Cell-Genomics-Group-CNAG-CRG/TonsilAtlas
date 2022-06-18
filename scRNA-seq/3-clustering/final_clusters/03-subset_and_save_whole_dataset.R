# This script subsets the level 1 Seurat object with the final
# cells that will compose the atlas and reprocesses it


# Load packages
library(Seurat)
library(tidyverse)
library(harmony)
library(here)


# Source functions
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Load data
seurat <- readRDS(path_to_tonsil_level_1)
final_cells_df <- readRDS(path_to_save_df_multi_myeloid)


# Subset cells
seurat <- subset(seurat, cells = final_cells_df$barcode)
seurat@meta.data <- final_cells_df


# Save
saveRDS(seurat, path_to_save_tonsil)
