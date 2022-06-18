# This script createss the supplementary tables associated with PDC, preT, preB


# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(openxlsx)


# Source utilities
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure2.R"))
path_to_obj <- here("data/raw_data_figures/seurat_obj_pdc_preB_preT_fig1.rds")


# Read data
seurat <- readRDS(path_to_obj)


# Find markers ALL
seurat$annotation_20220215 <- factor(
  seurat$annotation_20220215,
  levels = c("PDC", "IFN1+ PDC", "preB", "preT")
)
all_levels <- levels(seurat$annotation_20220215)
Idents(seurat) <- "annotation_20220215"
markers_all <- FindAllMarkers(
  seurat,
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
  here("results/paper/tables/supplementary_table_PDC_preB_preT.xlsx")
)
