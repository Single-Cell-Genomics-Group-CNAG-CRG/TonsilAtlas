# This script merges the demultiplexed Seurat objects into a single one


# Load packages
print("Loading packages...")
library(Seurat)
library(ggpubr)
library(tidyverse)


# Define paths
print("Defining paths...")
path_to_demultiplexed <- here::here("scRNA-seq/results/R_objects/demultiplexed/")
path_to_project_metadata <- here::here("scRNA-seq/1-cellranger_mapping/data/tonsil_atlas_metadata.csv")
path_to_donor_metadata <- here::here("data/tonsil_atlas_donor_metadata.csv")
path_to_save_obj <- here::here("scRNA-seq/results/R_objects/seurat_hashed_merged.rds")


# Read demultiplexed Seurat objects
print("Reading Seurat objects...")
files_to_load <- list.files(path_to_demultiplexed)
gem_ids <- files_to_load %>%
  str_remove("seurat_") %>%
  str_remove("_demultiplexed.rds")
files_to_load <- str_c(path_to_demultiplexed, files_to_load, sep = "")
tonsil_list <- purrr::map(files_to_load, readRDS)
names(tonsil_list) <- gem_ids


# Read project and donor metadata
print("Reading metadata...")
metadata_df <- read_csv(path_to_project_metadata)
donor_id_df <- read_csv(path_to_donor_metadata)


# Remove HTO assays and add some metadata
print("Pre-processing...")
metadata_df <- metadata_df[metadata_df$type == "hashed_cdna", ]
library_names <- metadata_df$library_name
donor_ids <- metadata_df$donor_id
names(library_names) <- metadata_df$gem_id
names(donor_ids) <- metadata_df$gem_id
donor_id_df <- as.data.frame(donor_id_df)
rownames(donor_id_df) <- donor_id_df$donor_id
tonsil_list <- purrr::map(gem_ids, function(x) {
  print(x)
  seurat_obj <- tonsil_list[[x]]
  seurat_obj[["HTO"]] <- NULL
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


# Merge Seurat objects
print("Merging Seurat objects...")
tonsil_hashed <- tonsil_list[[1]]
for (i in 2:length(tonsil_list)) {
  print(gem_ids[i])
  tonsil_hashed <- merge(x = tonsil_hashed, y = tonsil_list[[i]])
  tonsil_list[[i]] <- NA
}
rm(tonsil_list)


# Save
print("Saving...")
saveRDS(tonsil_hashed, path_to_save_obj)
