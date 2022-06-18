# This script finds markers specific for the epithelial populations


# Load packages
library(tidyverse)
library(Seurat)


# Load data
seurat <- readRDS(here::here("scRNA-seq/results/tonsil_atlas_all_cells_20210930.rds"))


# Subset to 3P
seurat <- subset(seurat, assay == "3P")


# Markers epithelial
epithelial_clusters <- unique(seurat$annotation_20210801[seurat$annotation_level_1 == "epithelial"])
Idents(seurat) <- "annotation_20210801"
markers_epithelial <- purrr::map(epithelial_clusters, function(x) {
  df <- FindMarkers(
    seurat,
    ident.1 = x,
    only.pos = TRUE,
    logfc.threshold = 1
  )
  df <- rownames_to_column(df, "gene")
  df
})
names(markers_epithelial) <- epithelial_clusters


# Save
saveRDS(
  markers_epithelial,
  here::here("scRNA-seq/results/R_objects/markers_final_cells/markers_epithelial_vs_whole_tonsil.rds")
)
