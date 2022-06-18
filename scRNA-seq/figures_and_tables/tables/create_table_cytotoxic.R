# This script plots the supplementary tables associated with cytotoxic cells


# Load packages
library(Seurat)
library(tidyverse)
library(here)
library(openxlsx)


# Source utilities
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure3.R"))


# Read data
cd8 <- readRDS(path_to_save_cd8)
ilc_nk <- readRDS(path_to_save_ilc_nk)


# Find markers ALL
ilc_nk$UMAP_1_20220215 <- ilc_nk$UMAP_1_20220215 + 10
ilc_nk$UMAP_2_20220215 <- ilc_nk$UMAP_2_20220215 + 10
levels(cd8$annotation_20220215)[levels(cd8$annotation_20220215) == "CXCR6+ RM CD8 T"] <- "RM CD8 activated T"
seurat_rna <- merge(x = cd8, y = ilc_nk)
seurat_rna$annotation_paper <- factor(
  seurat_rna$annotation_20220215,
  levels = names(colors_rna)
)
Idents(seurat_rna) <- "annotation_paper"
markers_all <- FindAllMarkers(
  seurat_rna,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
markers_all <- filter(markers_all, p_val_adj < 0.001)
markers_all_l <- map(names(colors_rna), function(x) {
  df <- markers_all[markers_all$cluster == x, ]
  df
})
names(markers_all_l) <- str_replace(names(colors_rna), ":", "/")


# Save
openxlsx::write.xlsx(
  markers_all_l,
  here("results/paper/tables/supplementary_table_cytotoxic_all.xlsx")
)
