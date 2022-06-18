# This script plots the supplementary tables associated with CD4 T cells


# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(openxlsx)


# Source utilities
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure2.R"))


# Read data
seurat_rna <- readRDS(path_to_save_cd4)
seurat_th <- readRDS(path_to_save_th)
path_to_atac_treg <- here("scATAC-seq/results/R_objects/Targetted_analysis/CD4_T/20220412_seurat_treg_atac.rds")
seurat_treg_atac <- readRDS(path_to_atac_treg)


# Rename
levels(seurat_rna$annotation_20220215)[levels(seurat_rna$annotation_20220215) == "GC-Tf-regs"] <- "Tfr"
levels(seurat_rna$annotation_20220215)[levels(seurat_rna$annotation_20220215) == "non-GC-Tf-regs"] <- "Eff-Tregs-IL32"
Idents(seurat_rna) <- "annotation_20220215"
levels(seurat_treg_atac$annotation_20220215)[levels(seurat_treg_atac$annotation_20220215) == "GC-Tf-regs"] <- "Tfr"
levels(seurat_treg_atac$annotation_20220215)[levels(seurat_treg_atac$annotation_20220215) == "non-GC-Tf-regs"] <- "Eff-Tregs-IL32"
Idents(seurat_treg_atac) <- "annotation_20220215"

# Find markers ALL
all_levels <- levels(seurat_rna$annotation_20220215)
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


# Find markers Naive and CM T
naive_cm_levels <- c("Naive", "CM Pre-non-Tfh", "CM PreTfh")
seurat_naive_cm <- subset(seurat_rna, idents = naive_cm_levels)
markers_naive_cm <- FindAllMarkers(
  seurat_naive_cm,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
markers_naive_cm <- filter(markers_naive_cm, p_val_adj < 0.001)
markers_naive_cm_l <- map(naive_cm_levels, function(x) {
  df <- markers_naive_cm[markers_naive_cm$cluster == x, ]
  df
})
names(markers_naive_cm_l) <- str_replace(naive_cm_levels, ":", "/")


# Find markers Th
th_levels <- levels(seurat_th$annotation_20220215)
markers_th <- FindAllMarkers(
  seurat_th,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
markers_th <- filter(markers_th, p_val_adj < 0.001)
markers_th_l <- map(th_levels, function(x) {
  df <- markers_th[markers_th$cluster == x, ]
  df
})
names(markers_th_l) <- str_replace(th_levels, ":", "/")

# Find markers Treg
treg_levels <- c("Eff-Tregs", "Eff-Tregs-IL32", "Tfr")
seurat_treg <- subset(seurat_rna, idents = treg_levels)
markers_treg <- FindAllMarkers(
  seurat_treg,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
markers_treg <- filter(markers_treg, p_val_adj < 0.001)
markers_treg_l <- map(treg_levels, function(x) {
  df <- markers_treg[markers_treg$cluster == x, ]
  df
})
treg_levels <- str_replace(treg_levels, ":", "/")
names(markers_treg_l) <- str_c("markers", treg_levels, sep = "_")


# Find accessible motifs Treg
markers_treg_atac <- FindAllMarkers(
  object = seurat_treg_atac,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
markers_treg_atac$motif_name <- unlist(
  seurat_treg_atac@assays$peaks_redefined@motifs@motif.names[rownames(markers_treg_atac)]
)
markers_treg_atac$motif <- rownames(markers_treg_atac)
markers_treg_atac <- filter(markers_treg_atac, p_val_adj < 0.001)
markers_treg_atac_l <- map(treg_levels, function(x) {
  df <- markers_treg_atac[markers_treg_atac$cluster == x, ]
  df
})
names(markers_treg_atac_l) <- str_c("motifs", treg_levels, sep = "_")


# Save
openxlsx::write.xlsx(
  markers_all_l,
  here("results/paper/tables/supplementary_table_CD4_T_all.xlsx")
)
openxlsx::write.xlsx(
  markers_naive_cm_l,
  here("results/paper/tables/supplementary_table_CD4_T_naive_cm.xlsx")
)
openxlsx::write.xlsx(
  markers_th_l,
  here("results/paper/tables/supplementary_table_CD4_T_Th.xlsx")
)
openxlsx::write.xlsx(
  c(markers_treg_l, markers_treg_atac_l),
  here("results/paper/tables/supplementary_table_CD4_T_Treg.xlsx")
)
