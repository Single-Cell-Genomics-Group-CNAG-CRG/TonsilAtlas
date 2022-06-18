# This script reads all non-hashed expression matrices and metges them into
# a single Seurat object


# Load packages
print("Loading packages...")
library(Seurat)
library(tidyverse)


# Define paths
print("Defining paths...")
path_to_data <- here::here("scRNA-seq/1-cellranger_mapping/projects")
path_to_project_metadata <- here::here("scRNA-seq/1-cellranger_mapping/data/tonsil_atlas_metadata.csv")
path_to_donor_metadata <- here::here("data/tonsil_atlas_donor_metadata.csv")
path_to_scrublet <- here::here("scRNA-seq/results/tables/scrublet")
path_to_save_obj <- here::here("scRNA-seq/results/R_objects/seurat_non_hashed_merged.rds")


# Read metadata and generate vector with paths to non-hashed matrices
print("Reading metadata...")
donor_id_df <- read_csv(path_to_donor_metadata)
metadata_df <- read_csv(path_to_project_metadata)
metadata_df <- filter(metadata_df, type == "not_hashed")
paths_to_data <- purrr::map2_chr(
  metadata_df$subproject,
  metadata_df$gem_id,
  function(subproject, gem_id) {
    path <- str_c(path_to_data, subproject, "jobs", gem_id, gem_id, "outs",
                  "filtered_feature_bc_matrix", sep = "/")
    path
  }
)
gem_ids <- metadata_df$gem_id
names(paths_to_data) <- gem_ids


# Read expression matrices and create Seurat objects
matrices <- purrr::map(paths_to_data, Read10X)
tonsil_list <- purrr::map(matrices, CreateSeuratObject)


# Add metadata
print("Pre-processing...")
library_names <- metadata_df$library_name
donor_ids <- metadata_df$donor_id
names(library_names) <- metadata_df$gem_id
names(donor_ids) <- metadata_df$gem_id
donor_id_df <- as.data.frame(donor_id_df)
rownames(donor_id_df) <- donor_id_df$donor_id
tonsil_list <- purrr::map(gem_ids, function(x) {
  print(x)
  seurat_obj <- tonsil_list[[x]]
  new_barcodes <- str_c(x, colnames(seurat_obj), sep = "_")
  seurat_obj <- RenameCells(seurat_obj, new.names = new_barcodes)
  seurat_obj$gem_id <- x
  seurat_obj$library_name <- library_names[x]
  seurat_obj$donor_id <- donor_ids[x]
  seurat_obj$hospital <- donor_id_df[seurat_obj$donor_id, "hospital"]
  seurat_obj$sex <- donor_id_df[seurat_obj$donor_id, "sex"]
  seurat_obj$age <- donor_id_df[seurat_obj$donor_id, "age"]
  seurat_obj$age_group <- donor_id_df[seurat_obj$donor_id, "age_group"]
  seurat_obj
})
names(tonsil_list) <- gem_ids


# Add scrublet doublet scores
list_scrublet_files <- list.files(path_to_scrublet)
paths_to_scrublet <- str_c(path_to_scrublet, list_scrublet_files, sep = "/")
tonsil_list <- purrr::map(gem_ids, function(x) {
  print(x)
  seurat_obj <- tonsil_list[[x]]
  scrublet_df <- paths_to_scrublet %>%
    str_subset(x) %>%
    read_csv(col_names = TRUE)
  scrublet_df$barcodes <- str_c(x, scrublet_df$barcodes, sep = "_")
  
  if (all(scrublet_df$barcodes == colnames(seurat_obj))) {
    warning("barcodes are equal")
    seurat_obj$scrublet_doublet_scores <- scrublet_df$scrublet_doublet_scores
    seurat_obj$scrublet_doublet_scores_scaled <- scale(
      scrublet_df$scrublet_doublet_scores,
      center = TRUE,
      scale = TRUE
    )
    seurat_obj$scrublet_predicted_doublet <- scrublet_df$scrublet_predicted_doublet
    
  } else{
    warning("barcodes are not equal")
  }
  
  seurat_obj
})
names(tonsil_list) <- gem_ids


# Merge Seurat objects
print("Merging Seurat objects...")
tonsil_non_hashed <- tonsil_list[[1]]
for (i in 2:length(tonsil_list)) {
  print(gem_ids[i])
  tonsil_non_hashed <- merge(x = tonsil_non_hashed, y = tonsil_list[[i]])
  tonsil_list[[i]] <- NA
}
rm(tonsil_list)


# Save
print("Saving...")
saveRDS(tonsil_non_hashed, path_to_save_obj)
