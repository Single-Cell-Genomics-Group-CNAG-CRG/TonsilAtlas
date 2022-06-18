# This script performs feature selection, dimensionality reduction and
# batch effect correction for a specific cell type (level 3)


# Load packages
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(tidyverse)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]


# Define parameters
path_to_level_2 <- here::here("scRNA-seq/results/R_objects/level_2/")
path_to_obj <- str_c(
  path_to_level_2,
  cell_type,
  "/",
  cell_type,
  "_clustered_filtered_level_2.rds",
  sep = ""
)
path_to_utils <- here::here("scRNA-seq/bin/utils.R")
path_to_level_3 <- here::here("scRNA-seq/results/R_objects/level_3/")
path_to_level_3_cell_type <- str_c(path_to_level_3, cell_type, sep = "")
path_to_save <- str_c(
  path_to_level_3_cell_type,
  "/",
  cell_type,
  "_integrated_level_3.rds",
  sep = ""
)


# Create directories (if they do not exist)
dir.create(path_to_level_3, showWarnings = FALSE)
dir.create(path_to_level_3_cell_type, showWarnings = FALSE)


# Source functions
source(path_to_utils)


# Load data
seurat <- readRDS(path_to_obj)


# Find consensus features
seurat$annotation_level_2 <- seurat$seurat_clusters
under_represented <- c("myeloid", "FDC", "PDC", "epithelial")
if (!(cell_type %in% under_represented)) {
  if (cell_type == "PC") {
    seurat <- subset(seurat, subset = assay != "5P")
  }
  seurat_list <- SplitObject(seurat, split.by = "assay")
  seurat_list <- seurat_list[c("3P", "multiome")]
  seurat_list <- purrr::map(
    seurat_list,
    FindVariableFeatures,
    nfeatures = 5000
  )
  hvg <- purrr::map(seurat_list, VariableFeatures)
  shared_hvg <- intersect(hvg$`3P`, hvg$multiome)
  print(length(shared_hvg))
} else {
  seurat <- subset(seurat, subset = assay == "3P")
}


# Dimensionality reduction and batch correction
if (!(cell_type %in% under_represented)) {
  seurat <- seurat %>%
      ScaleData(features = shared_hvg) %>%
      RunPCA(features = shared_hvg) %>%
      RunHarmony(group.by.vars = "assay", reduction = "pca", dims = 1:30)
} else {
  seurat <- seurat %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = "donor_id", reduction = "pca", dims = 1:30)
}


if (cell_type != "epithelial") {
  seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:30)
  seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:30)
} else {
  seurat <- RunUMAP(
    seurat,
    reduction = "harmony",
    dims = 1:20,
    n.neighbors = 10
  )
  seurat <- FindNeighbors(
    seurat,
    reduction = "harmony",
    dims = 1:20,
    k.param = 10
  )
}

  
# Save
saveRDS(seurat, path_to_save)