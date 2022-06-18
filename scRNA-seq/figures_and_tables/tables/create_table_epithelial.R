# This script creates all supplementary tables associated with epithelial cells

# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(openxlsx)


# Source utilities
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure4.R"))


# Read data
seurat_rna <- readRDS(path_to_save_epithelial)


# Markers all
all_levels <- names(colors_epi)
seurat_rna$annotation_20220215 <- factor(
  seurat_rna$annotation_20220215,
  levels = all_levels
)
Idents(seurat_rna) <- "annotation_20220215"
markers_all <- FindAllMarkers(
  seurat_rna,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
markers_all <- filter(markers_all, p_val_adj < 0.001)
markers_all_l <- map(all_levels, function(x) {
  df <- markers_all[markers_all$cluster == x, ]
  df
})
names(markers_all_l) <- str_replace(all_levels, ":", "/")


# Save
openxlsx::write.xlsx(
  markers_all_l,
  here("results/paper/tables/supplementary_table_epithelial_markers_all.xlsx")
)
