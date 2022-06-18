# This script includes important metadata in the Seurat object and filters out
# doublets that were predicted by Scrublet in both ATAC and RNA.


# Load packages
library(Seurat)
library(Signac)
library(tidyverse)


# Define paths
path_to_data <- here::here("multiome/results/R_objects/3.tonsil_multiome_filtered_combined.rds")
path_to_multiome_metadata <- here::here("multiome/1-cellranger_mapping/data/tonsil_atlas_metadata_multiome.csv")
path_to_donor_metadata <- here::here("data/tonsil_atlas_donor_metadata.csv")
path_to_scrublet_dfs <- here::here("multiome/results/tables/scrublet/")
path_to_save <- here::here("multiome/results/R_objects/4.tonsil_multiome_filtered_combined_with_metadata.rds")


# Load data
tonsil_multiome <- readRDS(path_to_data)
multiome_metadata <- read_csv(path_to_multiome_metadata)
donor_metadata <- read_csv(path_to_donor_metadata)
scrublet_files <- list.files(path_to_scrublet_dfs)
scrublet_files <- str_c(path_to_scrublet_dfs, scrublet_files, sep = "")
all_scrublet <- purrr::map(scrublet_files, read_csv)
names(all_scrublet) <- list.files(path_to_scrublet_dfs) %>% 
  str_remove(".csv") %>% 
  str_remove("scrublet_doublet_prediction_") %>%
  str_remove("sparse_matrix_")
new_colnames_atac <- c("scrublet_doublet_scores_atac", "scrublet_predicted_doublet_atac")
colnames(all_scrublet$atac_with_BCLL_2)[2:3] <- new_colnames_atac
colnames(all_scrublet$atac_without_BCLL_2)[2:3] <- new_colnames_atac


# Include metadata
multiome_metadata_sub <- multiome_metadata %>%
  group_by(gem_id) %>%
  dplyr::filter(row_number(`gem_id`) == 1)
donor_ids <- multiome_metadata_sub$donor_id
names(donor_ids) <- multiome_metadata_sub$gem_id
tonsil_multiome@meta.data$donor_id <- donor_ids[tonsil_multiome$gem_id]
sex_vec <- donor_metadata$sex
age_vec <- donor_metadata$age
age_group_vec <- donor_metadata$age_group
hospital_vec <- donor_metadata$hospital
names(sex_vec) <- donor_metadata$donor_id
names(age_vec) <- donor_metadata$donor_id
names(age_group_vec) <- donor_metadata$donor_id
names(hospital_vec) <- donor_metadata$donor_id
tonsil_multiome@meta.data$sex <- sex_vec[tonsil_multiome$donor_id]
tonsil_multiome@meta.data$age <- age_vec[tonsil_multiome$donor_id]
tonsil_multiome@meta.data$age_group <- age_group_vec[tonsil_multiome$donor_id]
tonsil_multiome@meta.data$hospital <- hospital_vec[tonsil_multiome$donor_id]
tonsil_multiome@meta.data$assay <- "multiome"


# Include doublet scores
bcll_2_scrublet_df <- left_join(
  x = all_scrublet$rna_with_BCLL_2,
  y = all_scrublet$atac_with_BCLL_2,
  by = "barcodes"
)
bcll_2_scrublet_df <- bcll_2_scrublet_df %>%
  mutate(
    scrublet_doublet_scores_scaled = scale(scrublet_doublet_scores),
    scrublet_doublet_scores_scaled_atac = scale(scrublet_doublet_scores_atac)
  )
non_bcll_2_scrublet_df <- left_join(
  x = all_scrublet$rna_without_BCLL_2,
  y = all_scrublet$atac_without_BCLL_2,
  by = "barcodes"
)
non_bcll_2_scrublet_df <- non_bcll_2_scrublet_df %>%
  mutate(
    scrublet_doublet_scores_scaled = scale(scrublet_doublet_scores),
    scrublet_doublet_scores_scaled_atac = scale(scrublet_doublet_scores_atac)
  )
scrublet_df <- bind_rows(list(bcll_2_scrublet_df, non_bcll_2_scrublet_df))
scrublet_df <- as.data.frame(scrublet_df)
rownames(scrublet_df) <- scrublet_df$barcodes
scrublet_df <- scrublet_df[colnames(tonsil_multiome), ]
tonsil_multiome@meta.data$scrublet_doublet_scores <- scrublet_df$scrublet_doublet_scores
tonsil_multiome@meta.data$scrublet_predicted_doublet <- scrublet_df$scrublet_predicted_doublet
tonsil_multiome@meta.data$scrublet_doublet_scores_scaled <- scrublet_df$scrublet_doublet_scores_scaled
tonsil_multiome@meta.data$scrublet_doublet_scores_atac <- scrublet_df$scrublet_doublet_scores_atac
tonsil_multiome@meta.data$scrublet_predicted_doublet_atac <- scrublet_df$scrublet_predicted_doublet_atac
tonsil_multiome@meta.data$scrublet_doublet_scores_scaled_atac <- scrublet_df$scrublet_doublet_scores_scaled_atac


# Save
saveRDS(tonsil_multiome, path_to_save)
