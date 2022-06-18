# This script merges and save the Seurat objects of dendritic cells (DC) and naive CD4 T cells.
# We will use them as a positive control for the cell-cell interaction analysis. That is, if we
# recover the interactions that have been clearly identified in the literature, then we will be
# in a good position to scale the analysis to more cell types.


# Load packages
library(Matrix)
library(Seurat)
library(tidyverse)


# Parameters
path_to_cd4 <- here::here("scRNA-seq/results/R_objects/level_4/CD4_T/CD4_T_integrated_level_4.rds")
path_to_myeloid <- here::here("scRNA-seq/results/R_objects/level_3/myeloid/myeloid_clustered_level_3_with_pre_freeze.rds")
path_to_save_barcodes <- here::here("")
path_to_save_features <- here::here("")
path_to_save_counts <- here::here(".mtx")
path_to_save_annotation <- here::here("")


# Load functions
source(here::here("scRNA-seq/bin/utils.R"))


# Read data
cd4 <- readRDS(path_to_cd4)
myeloid <- readRDS(path_to_myeloid)


# Subset
naive_cd4 <- subset(cd4, idents = "Naive")
Idents(myeloid) <- "annotation_pre_freeze"
dc <- subset(myeloid, idents = c("DC1", "DC2", "aDC", "aDC2"))


# Merge
dc$annotation_level_3 <- dc$annotation_pre_freeze
seurat <- merge(x = naive_cd4, y = dc)
seurat <- subset(seurat, subset = assay == "3P")


# Save
counts <- seurat[["RNA"]]@counts
features <- rownames(seurat)
barcodes <- colnames(seurat)
metadata <- seurat@meta.data[, "annotation_level_3", drop = FALSE]
metadata <- rownames_to_column(metadata, var = "barcode")
write_csv(metadata, file = path_to_save_annotation)
writeMM(t(counts), path_to_save_counts) # Remember scanpy works with transposed matrices
write(barcodes, path_to_save_barcodes)
write(features, path_to_save_features)
