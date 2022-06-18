# This script computes the Local Inverse Simpson Index (LISI) to quantify the
# effect of different confounders in the tonsil atlas pre- and post-batch
# effect correction


# Load packages
library(tidyverse)
library(lisi)


# Define paths
path_tmp_dir <- here::here("scRNA-seq/2-QC/5-batch_effect_correction/3-data_integration_multiome/tmp/")
path_to_unintegrated_dimred <- str_c(path_tmp_dir, "batch_uncorrected_pca.rds", sep = "")
path_to_integrated_dimred <- str_c(path_tmp_dir, "batch_corrected_pca.rds", sep = "")
path_to_confounders_df <- str_c(path_tmp_dir, "confounders_df.rds", sep = "") 
path_to_save_lisi_scores <- str_c(path_tmp_dir, "lisi_scores.rds", sep = "") 


# Load data
pca_unintegrated <- readRDS(path_to_unintegrated_dimred)
pca_integrated <- readRDS(path_to_integrated_dimred)
confounders_df <- readRDS(path_to_confounders_df)


# Compute LISI
dim_red_mats <- list(pca_unintegrated, pca_integrated)
names(dim_red_mats) <- c("unintegrated", "integrated")
confounders_df$is_hashed[is.na(confounders_df$is_hashed)] <- "NA"
lisi_scores <- purrr::map(dim_red_mats, function(mat) {
  scores <- compute_lisi(
    X = mat,
    meta_data = confounders_df,
    label_colnames = colnames(confounders_df)
  )
  scores
})
lisi_scores <- bind_rows(lisi_scores, .id = "is_integrated")


# Save
saveRDS(lisi_scores, path_to_save_lisi_scores)

