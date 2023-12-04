# This script splits the merged tonsil atlas (discovery + validation) into
# one Seurat object per compartment, and saves them independently


# Load packages
library(Seurat)
library(harmony)
library(tidyverse)
library(here)
library(glue)


# Load functions
source(here("scRNA-seq/bin/SLOcatoR_functions.R"))


# Read data
seurat <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/3.1-seurat_object_discovery_validation_cohorts_integrated.rds"))
metadata_discovery <- readRDS(
  here("scRNA-seq/results/R_objects/seurat_objects_revision/cell_metadata_discovery_cohort.rds")
)
metadata_validation <- readRDS(
  here("scRNA-seq/results/R_objects/seurat_objects_revision/cell_metadata_validation_cohort.rds")
)
donor_metadata <- read_csv(here("data/tonsil_atlas_donor_metadata.csv"))


# Harmonize metadata
seurat$library_name[is.na(seurat$library_name)] <- seurat$gem_id[is.na(seurat$library_name)] 
seurat$cohort_type <- seurat$type
seurat$cohort_type <- str_remove(seurat$cohort_type, "(_scRNA|_multiome)")
donor_metadata <- donor_metadata[, c("donor_id", "cause_for_tonsillectomy")]
metadata_all <- left_join(seurat@meta.data, donor_metadata, by = "donor_id")
if (all(colnames(seurat) == metadata_all$barcode)) {
  seurat$cause_for_tonsillectomy <- metadata_all$cause_for_tonsillectomy
}
seurat@meta.data$batch <- NULL
seurat@meta.data$type2 <- NULL
colnames(seurat@meta.data)[colnames(seurat@meta.data) == "label"] <- "annotation_level_1"
colnames(seurat@meta.data)[colnames(seurat@meta.data) == "annotation_probability"] <- "annotation_level_1_probability"
colnames(seurat@meta.data)[colnames(seurat@meta.data) == "UMAP1"] <- "UMAP_1_level_1"
colnames(seurat@meta.data)[colnames(seurat@meta.data) == "UMAP2"] <- "UMAP_2_level_1"
metadata_all2 <- left_join(seurat@meta.data, metadata_discovery, "barcode")
if (all(colnames(seurat) == metadata_all2$barcode)) {
  seurat$is_hashed <- metadata_all2$is_hashed
  seurat$is_hashed[seurat$cohort_type == "validation"] <- "not_hashed"
  seurat$pDNN_hashing <- metadata_all2$pDNN_hashing
  seurat$pDNN_scrublet <- metadata_all2$pDNN_scrublet
  seurat$pDNN_union <- metadata_all2$pDNN_union
  seurat$scrublet_doublet_scores <- metadata_all2$scrublet_doublet_scores
  seurat$scrublet_predicted_doublet <- metadata_all2$scrublet_predicted_doublet
  seurat$annotation_20220215 <- metadata_all2$annotation_20220215
  seurat$annotation_figure_1 <- metadata_all2$annotation_figure_1
  seurat$UMAP_1_20220215 <- metadata_all2$UMAP_1_20220215
  seurat$UMAP_2_20220215 <- metadata_all2$UMAP_2_20220215
}
seurat$doublet_score_scDblFinder <- NA
seurat$doublet_score_scDblFinder[rownames(metadata_validation)] <- metadata_validation$doublet_score
seurat$preservation <- ifelse(seurat$assay == "3P", "fresh", "frozen")
seurat$preservation[seurat$donor_id %in% c("BCLL-25-T", "BCLL-26-T")] <- "frozen"
seurat$preservation[str_detect(seurat$gem_id, "frozen")] <- "frozen"


# Recompute cell cycle scores
seurat <- CellCycleScoring(
  seurat,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE
)


# Split
seurat_list <- SplitObject(seurat, split.by = "annotation_level_1")


# Save
for (i in 1:length(seurat_list)) {
  x <- names(seurat_list)[i]
  print(x)
  if (x != "preTC") {
    shared_hvg <- find_assay_specific_features(seurat_list[[x]], assay_var = "type")
    print(glue("The number of shared HVG is: {length(shared_hvg)}"))
    seurat_list[[x]] <- integrate_assays(
      seurat_list[[x]],
      assay_var = "type",
      shared_hvg = shared_hvg
    )
    seurat_list[[x]] <- RunUMAP(seurat_list[[x]], dims = 1:25, reduction = "harmony")
  }
  seurat_list[[x]]@meta.data$type <- NULL
  saveRDS(
    seurat_list[[x]],
    here(glue("scRNA-seq/results/R_objects/seurat_objects_revision/4.{i}-{x}_seurat_object_discovery_validation_cohorts.rds"))
  )
}

