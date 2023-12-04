# This script merges the different Seurat objects we are aiming to add for
# the revision


# Load packages
print("Loading packages...")
library(Seurat)
library(Signac)
library(tidyverse)
library(here)


# Read data
print("Reading data...")
seurat_newcastle <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_newcastle/1-seurat_object_newcastle_filtered.rds"))
seurat_bcn_rna <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/1-seurat_object_barcelona_filtered.rds"))
seurat_bcn_multiome <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/1-seurat_object_barcelona_multiome_filtered.rds"))
donor_metadata <- read_csv(here("data/tonsil_atlas_donor_metadata.csv"))
sequencing_metadata_rna <- read_csv(here("scRNA-seq/1-cellranger_mapping/data/tonsil_atlas_metadata.csv"))
sequencing_metadata_multiome <- read_csv(here("multiome/1-cellranger_mapping/data/tonsil_atlas_metadata_multiome.csv"))


# Merge multiome
print("Merging multiome...")
n_cells_list <- map(seurat_bcn_multiome, \(seurat_obj) {
  n_cells <- Matrix::rowSums(seurat_obj[["RNA"]]@counts > 0)
  n_cells
})
seurat_bcn_multiome <- map2(seurat_bcn_multiome, n_cells_list, \(seurat_obj, n_cells) {
  selected_genes <- names(n_cells)[n_cells > 5]
  seurat_obj_atac <- seurat_obj[["ATAC"]]
  seurat_obj <- subset(seurat_obj, features = selected_genes)
  seurat_obj[["ATAC"]] <- seurat_obj_atac
  seurat_obj
})
merged_multiome <- Reduce(merge, seurat_bcn_multiome)


# Metadata
print("Including metadata...")
sequencing_metadata_multiome <- sequencing_metadata_multiome[sequencing_metadata_multiome$gem_id %in% unique(merged_multiome$gem_id), ]
sequencing_metadata_multiome <- sequencing_metadata_multiome[sequencing_metadata_multiome$type == "RNA", ]
new_metadata_multiome <- left_join(
  merged_multiome@meta.data,
  sequencing_metadata_multiome,
  by = "gem_id"
)
new_metadata_multiome <- left_join(
  new_metadata_multiome,
  donor_metadata,
  by = "donor_id"
)
new_metadata_multiome$barcode <- colnames(merged_multiome)
new_metadata_multiome <- as.data.frame(new_metadata_multiome)
rownames(new_metadata_multiome) <- new_metadata_multiome$barcode
merged_multiome@meta.data <- new_metadata_multiome
saveRDS(merged_multiome, here("scRNA-seq/results/R_objects/seurat_objects_revision/2.2-seurat_object_barcelona_multiome_merged.rds"))


# Merge all
print("Merging all..")
merged_multiome[["ATAC"]] <- NULL
merged <- merge(x = seurat_bcn_rna, y = list(seurat_newcastle, merged_multiome))


# Save
print("Saving...")
saveRDS(merged, here("scRNA-seq/results/R_objects/seurat_objects_revision/2.1-seurat_object_barcelona_newcastle_multiome_merged.rds"))