# This script merges the hashed and non-hashed Seurat objects into a single one


# Load packages
print("Loading packages...")
library(Seurat)
library(tidyverse)


# Define paths
print("Defining paths...")
path_to_hashed <- here::here("scRNA-seq/results/R_objects/seurat_hashed_merged_filtered.rds")
path_to_not_hashed <- here::here("scRNA-seq/results/R_objects/seurat_non_hashed_merged_filtered.rds")
path_to_save_obj <- here::here("scRNA-seq/results/R_objects/seurat_merged_all.rds")


# Read Seurat objects
print("Reading Seurat objects...")
tonsil_hashed <- readRDS(path_to_hashed)
tonsil_not_hashed <- readRDS(path_to_not_hashed)


# Homogenize metadata
print("Homogenizing metadata...")
selected_variables <- c("nCount_RNA", "nFeature_RNA", "pct_mt", "pct_ribosomal",
                        "gem_id", "library_name", "donor_id", "sex", "age",
                        "age_group", "hospital", "scrublet_doublet_scores",
                        "scrublet_doublet_scores_scaled", "scrublet_predicted_doublet",
                        "HTO_classification.global", "has_high_lib_size")
print("Metadata variables tonsil_hashed")
print(colnames(tonsil_hashed@meta.data))
tonsil_hashed@meta.data <- tonsil_hashed@meta.data[, selected_variables]
tonsil_not_hashed$HTO_classification.global <- "NA"
print("Metadata variables tonsil_not_hashed:")
print(colnames(tonsil_not_hashed@meta.data))
tonsil_not_hashed@meta.data <- tonsil_not_hashed@meta.data[, selected_variables]
tonsil_hashed$is_hashed <- "hashed"
tonsil_not_hashed$is_hashed <- "not_hashed"


# Merge
print("Merging...")
tonsil <- merge(x = tonsil_hashed, y = tonsil_not_hashed)
tonsil


# Save
print("Saving...")
saveRDS(tonsil, path_to_save_obj)

