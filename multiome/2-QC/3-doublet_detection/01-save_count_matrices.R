# This script loads the Seurat object derived from the multiome data,
# fetches the counts matrices (ATAC + RNA), divides it into BCLL-2 and non-BCLL-2
# and saves them to effectively run scrublet on them


# Load packages
library(Signac)
library(Seurat)
library(Matrix)
library(tidyverse)


# Define paths
path_to_data <- here::here("multiome/results/R_objects/3.tonsil_multiome_filtered_combined.rds")
path_to_tmp_dir <- here::here("multiome/2-QC/3-doublet_detection/tmp/")


# Load data
tonsil_multiome <- readRDS(path_to_data)


# Create matrices
bcll_2_gem_ids <- c("kmbfo1ab_ie02b4ny", "ryh4el3i_biv0w7ca", "bs2e7lr7_mdfwypvz")
bcll_2_cells <- colnames(tonsil_multiome)[tonsil_multiome$gem_id %in% bcll_2_gem_ids]
non_bcll_2_cells <- colnames(tonsil_multiome)[!(tonsil_multiome$gem_id %in% bcll_2_gem_ids)]
rna_mat_bcll_2 <- tonsil_multiome[["RNA"]]@counts[, bcll_2_cells]
rna_mat_non_bcll_2 <- tonsil_multiome[["RNA"]]@counts[, non_bcll_2_cells]
atac_mat_bcll_2 <- tonsil_multiome[["peaks"]]@counts[, bcll_2_cells]
atac_mat_non_bcll_2 <- tonsil_multiome[["peaks"]]@counts[, non_bcll_2_cells]



# Save
dir.create(path_to_tmp_dir, showWarnings = FALSE)
mats_list <- list(
  rna_sparse_matrix_with_BCLL_2 = rna_mat_bcll_2,
  rna_sparse_matrix_without_BCLL_2 = rna_mat_non_bcll_2,
  atac_sparse_matrix_with_BCLL_2 = atac_mat_bcll_2,
  atac_sparse_matrix_without_BCLL_2 = atac_mat_non_bcll_2
)
for (x in names(mats_list)) {
  path_to_subdir <- str_c(path_to_tmp_dir, x, sep = "")
  dir.create(path_to_subdir, showWarnings = FALSE)
  path_save_mat <- str_c(path_to_subdir, "matrix.mtx", sep = "/")
  path_save_features <- str_c(path_to_subdir, "features.tsv", sep = "/")
  path_save_cell_barcodes <- str_c(path_to_subdir, "barcodes.tsv", sep = "/")
  writeMM(mats_list[[x]], path_save_mat)
  write(x = rownames(mats_list[[x]]), file = path_save_features)
  write(x = colnames(mats_list[[x]]), file = path_save_cell_barcodes)
}
