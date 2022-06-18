# This script load the integrated Seurat object and calculates de KNN graph


# Load packages
library(Seurat)


# Define paths
path_to_obj <- here::here("scRNA-seq/results/R_objects/tonsil_rna_integrated_with_multiome.rds")
path_to_save <- here::here("scRNA-seq/results/R_objects/tonsil_rna_integrated_knn.rds")


# Load data
tonsil <- readRDS(path_to_obj)


# Calculate the K-nearest neighbor graph
tonsil <- FindNeighbors(tonsil, reduction = "harmony", dims = 1:30)


# Save
saveRDS(tonsil, path_to_save)
