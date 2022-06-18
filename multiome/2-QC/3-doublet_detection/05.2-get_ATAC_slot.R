# This script loads the multiome Seurat object and saves only the RNA assay to
# merge it later with the scRNA-seq object


# Load packages
library(Seurat)
library(Signac)


# Define paths
path_to_data <- here::here("multiome/results/R_objects/4.tonsil_multiome_filtered_combined_with_metadata.rds")
path_to_save <- here::here("multiome/results/R_objects/5.tonsil_multiome_peaks_only.rds")


# Load data
tonsil_multiome <- readRDS(path_to_data)


# Exclude RNA
DefaultAssay(tonsil_multiome) <- "peaks"
tonsil_multiome[["RNA"]] <- NULL


# Save
saveRDS(tonsil_multiome, path_to_save)
