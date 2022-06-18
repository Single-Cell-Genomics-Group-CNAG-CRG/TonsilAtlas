# This script loads the annotated Seurat object (level 1) and creates a cell
# type-specific object and saves them individually.


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_obj <- here::here("scATAC-seq/results/R_objects/8.3.tonsil_peakcalling_annotation_level1_signature.rds")
path_to_save <- here::here("scATAC-seq/results/R_objects/level_2/")


# Load data
tonsil <- readRDS(path_to_obj)
tonsil


# Add current UMAP coords to metadata
umap_1 <- tonsil@reductions$umap@cell.embeddings[, "UMAP_1"]
umap_2 <- tonsil@reductions$umap@cell.embeddings[, "UMAP_2"]
tonsil$UMAP_1_level_1 <- umap_1
tonsil$UMAP_2_level_1 <- umap_2


# Split and save
tonsil_list <- SplitObject(tonsil, split.by = "annotation_level_1")
dir.create(path_to_save, showWarnings = FALSE)
for (cell_type in names(tonsil_list)) {
  print(cell_type)
  path_to_save_cell_type <- str_c(path_to_save, cell_type)
  dir.create(path_to_save_cell_type, showWarnings = FALSE)
  path_to_save_final <- str_c(
    path_to_save_cell_type,
    "/",
    cell_type,
    "_subsetted_level_2.rds",
    sep = ""
  )
  saveRDS(tonsil_list[[cell_type]], path_to_save_final)
}