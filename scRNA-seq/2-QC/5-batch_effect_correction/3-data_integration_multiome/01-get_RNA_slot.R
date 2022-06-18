# This script loads the multiome Seurat object and saves only the RNA assay to
# merge it later with the scRNA-seq object


# Load packages
library(Seurat)
library(Signac)


# Define paths
path_to_data <- here::here("multiome/results/R_objects/6.tonsil_multiome_integrated_using_wnn_no_doublets.rds")
path_to_save <- here::here("scRNA-seq/results/R_objects/7.tonsil_multiome_rna_only.rds")


# Load data
tonsil_multiome <- readRDS(path_to_data)


# Exclude ATAC
DefaultAssay(tonsil_multiome) <- "RNA"
tonsil_multiome[["peaks"]] <- NULL
tonsil_multiome@reductions[["pca"]] <- NULL
tonsil_multiome@reductions[["harmony_RNA"]] <- NULL
tonsil_multiome@reductions[["umap.atac"]] <- NULL
tonsil_multiome@reductions[["umap.rna"]] <- NULL
#tonsil_multiome@reductions[["wnn.umap"]] <- NULL


# Save
saveRDS(tonsil_multiome, path_to_save)
