# This script loads both the multiome and the scATAC-seq Seurat objects,
# harmonizes their metadata and merges them into a single object.

# Load packages
library(Seurat)
library(tidyverse)
library(plyr)
library(Signac)
library(GenomicRanges)
library(ggpubr)

## Functions
peaks_quantification <- function(seurat_filtered, peaks = combined.peaks){ 
  counts <- FeatureMatrix(
    fragments = Fragments(seurat_filtered),
    features = combined.peaks,
    cells = colnames(seurat_filtered))
    seurat_filtered[["peaks"]] <- CreateChromatinAssay(
    counts, 
    fragments = Fragments(seurat_filtered))
  return(seurat_filtered)}

remove_assays <- function(seurat_object){
  DefaultAssay(seurat_object) <- "peaks"
  seurat_object@assays[["ATAC"]] <- NULL
  return(seurat_object)
}

# Define paths
path_to_atac <- here::here("scATAC-seq/results/R_objects/4.tonsil_aggregated_harmony.rds")
path_to_multiome <- here::here("scATAC-seq/results/R_objects/5.tonsil_multiome_peaks_only_integrated_using_wnn_no_doublets.rds")
path_to_save <- here::here("scATAC-seq/results/R_objects/6.tonsil_atac_merged_with_multiome_missing_integration.rds")

# Load data
tonsil_atac <- readRDS(path_to_atac)
tonsil_multiome <- readRDS(path_to_multiome)

names(tonsil_atac@meta.data)[names(tonsil_atac@meta.data) == "scrublet_doublet_scores"] <- "scrublet_doublet_scores_atac"
names(tonsil_atac@meta.data)[names(tonsil_atac@meta.data) == "scrublet_predicted_doublet"] <- "scrublet_predicted_doublet_atac"

# Rename tonsil_multiome assay
tonsil_multiome <- RenameAssays(tonsil_multiome, peaks = "ATAC")

# Delete RNA-specific variables
tonsil_multiome$orig.ident <- NULL
tonsil_multiome$nCount_RNA <- NULL
tonsil_multiome$nFeature_RNA <- NULL
tonsil_multiome$pct_mt <- NULL
tonsil_multiome$pct_ribosomal <- NULL
tonsil_multiome$scrublet_predicted_doublet <- NULL
tonsil_multiome$scrublet_doublet_scores <- NULL
tonsil_multiome$scrublet_doublet_scores_scaled <- NULL
tonsil_multiome$scrublet_doublet_scores_scaled_atac <- NULL

# Create a unified set of peaks
combined.peaks <- UnifyPeaks(object.list = list(tonsil_atac,tonsil_multiome), mode = "reduce")   

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)

ggviolin(peakwidths,add = "boxplot",fill = "gray") + scale_y_log10() + 
  geom_hline(yintercept = c(20,10000), linetype='dashed', col = 'black')

combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
tonsil_unified_peaks <- lapply(list(tonsil_atac,tonsil_multiome), peaks_quantification)


# Removing ATAC assay to reduce the size of the final seurat object,
tonsil_unified_peaks <- lapply(tonsil_unified_peaks, remove_assays)

## Merging dataset
object1 <- tonsil_unified_peaks[[1]]
object1$assay <- "scATAC"
object2 <- tonsil_unified_peaks[[2]]
object2$assay <- "multiome"

Idents(object1) <- "gem_id"
Idents(object2) <- "gem_id"

tonsil <- merge(object1, y = object2)

# Save
saveRDS(tonsil, path_to_save)
