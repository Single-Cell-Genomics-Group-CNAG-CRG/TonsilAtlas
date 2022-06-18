# This script performs feature selection, dimensionality reduction and
# batch effect correction for a specific cell type (level 2)


# Load packages
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(tidyverse)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]


# Define parameters
path_to_dir <- here::here("scRNA-seq/results/R_objects/level_2/")
path_to_obj <- str_c(
  path_to_dir,
  cell_type,
  "/",
  cell_type,
  "_subsetted_level_2.rds",
  sep = ""
)
path_to_utils <- here::here("scRNA-seq/bin/utils.R")
path_to_save <- str_c(
  here::here("scRNA-seq/results/R_objects/level_2/"),
  cell_type,
  "/",
  cell_type,
  "_integrated_level_2.rds",
  sep = ""
)


# Source functions
source(path_to_utils)


# Load data
seurat <- readRDS(path_to_obj)


# Find consensus features
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


# Dimensionality reduction and batch correction
seurat <- seurat %>%
  ScaleData(features = shared_hvg) %>%
  RunPCA(features = shared_hvg) %>%
  RunHarmony(group.by.vars = "assay", reduction = "pca", dims = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)


# Save
saveRDS(seurat, path_to_save)


