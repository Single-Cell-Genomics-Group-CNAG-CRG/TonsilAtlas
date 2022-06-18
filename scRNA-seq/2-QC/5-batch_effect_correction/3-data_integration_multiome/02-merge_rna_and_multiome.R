# This script loads both the multiome and the scRNA-seq Seurat objects,
# harmonizes their metadata and merges them into a single object.


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_rna <- here::here("scRNA-seq/results/R_objects/seurat_merged_with_king_et_al_integrated.rds")
path_to_multiome <- here::here("scRNA-seq/results/R_objects/7.tonsil_multiome_rna_only.rds")
path_to_save <- here::here("scRNA-seq/results/R_objects/tonsil_rna_merged_with_multiome_missing_integration.rds")


# Load data
tonsil_rna <- readRDS(path_to_rna)
tonsil_multiome <- readRDS(path_to_multiome)


# Delete ATAC-specific variables
tonsil_multiome$orig.ident <- NULL
tonsil_multiome$nCount_ATAC <- NULL
tonsil_multiome$nFeature_ATAC <- NULL
tonsil_multiome$nCount_peaks <- NULL
tonsil_multiome$nFeature_peaks <- NULL
tonsil_multiome$nucleosome_signal <- NULL
tonsil_multiome$nucleosome_percentile <- NULL
tonsil_multiome$TSS.enrichment <- NULL
tonsil_multiome$TSS.percentile <- NULL
tonsil_multiome$high.tss <- NULL
tonsil_multiome$scrublet_predicted_doublet_atac <- NULL
tonsil_multiome$scrublet_doublet_scores_atac <- NULL
tonsil_multiome$scrublet_doublet_scores_scaled_atac <- NULL
tonsil_multiome$RNA_snn_res.1.25 <- NULL
tonsil_multiome$RNA_snn_res.1.5 <- NULL
tonsil_multiome$RNA_snn_res.1.75 <- NULL
tonsil_multiome$RNA_snn_res.2 <- NULL
tonsil_multiome$seurat_clusters <- NULL
tonsil_multiome$is_doublet <- NULL


# Calculate cell cycle scores (multiome)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
tonsil_multiome <- CellCycleScoring(
  tonsil_multiome,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)
tonsil_multiome$CC.Difference <- tonsil_multiome$S.Score - tonsil_multiome$G2M.Score


# Set IG columns to NA in our dataset
ig_columns <- str_subset(colnames(tonsil_rna@meta.data), "^IG")
for (x in ig_columns) {
  tonsil_multiome@meta.data[[x]] <- NA
}
tonsil_multiome$chain_status <- NA
tonsil_multiome$isotype <- NA
tonsil_multiome$sequence_id <- NA
tonsil_multiome$clone_id <- NA


# Set annotation columns to NA in our dataset
tonsil_multiome$cell_type <- "unannotated"
tonsil_multiome$lineage <- "unannotated"
tonsil_multiome$subset <- "unannotated"
tonsil_multiome$MBC_subset <- "unannotated"


# Set columns of tonsil_multiome to NA
tonsil_multiome$is_hashed <- NA
tonsil_multiome$HTO_classification.global <- NA
tonsil_multiome$has_high_lib_size <- NA
tonsil_multiome$pDNN_hashing <- NA
tonsil_multiome$pDNN_scrublet <- NA
tonsil_multiome$pDNN_union <- NA
tonsil_multiome$doublet_finder_predicted_doublet <- NA
tonsil_multiome$status <- NA


# Merge
Idents(tonsil_rna) <- "gem_id"
Idents(tonsil_multiome) <- "gem_id"
tonsil <- merge(x = tonsil_rna, y = tonsil_multiome)


# Save
saveRDS(tonsil, path_to_save)
