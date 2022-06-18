# This script creates all supplementary tables associated with myeloid cells



# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(openxlsx)


# Source utilities
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure4.R"))
alpha <- 0.05


# Read data
seurat_rna <- readRDS(path_to_save_myeloid)
gsea_results <- readRDS(here("scRNA-seq/results/R_objects/gsea_list_slancytes.rds"))


# Markers all
all_levels <- names(colors_myel)
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


# Markers slan+ cells
slan_levels <- c("IL7R MMP12 macrophages", "ITGAX ZEB2 macrophages",
                 "C1Q HLA macrophages", "SELENOP FUCA1 PTGDS macrophages")
seurat_slan <- subset(seurat_rna, idents = slan_levels)
markers_slan <- FindAllMarkers(
  seurat_slan,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
markers_slan <- filter(markers_slan, p_val_adj < 0.001)
markers_slan_l <- map(slan_levels, function(x) {
  df <- markers_slan[markers_slan$cluster == x, ]
  df
})
names(markers_slan_l) <- str_replace(slan_levels, ":", "/")


# Markers DC
dc_levels <- c("DC1 precursor", "DC1 mature", "DC2", "DC3", "DC4", "DC5",
               "IL7R DC")
seurat_dc <- subset(seurat_rna, idents = dc_levels)
markers_dc <- FindAllMarkers(
  seurat_dc,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
markers_dc <- filter(markers_dc, p_val_adj < 0.001)
markers_dc_l <- map(dc_levels, function(x) {
  df <- markers_dc[markers_dc$cluster == x, ]
  df
})
names(markers_dc_l) <- str_replace(dc_levels, ":", "/")


# table GO terms
gsea_sorted <- purrr::map(gsea_results, function(x) {
  df <- x@result %>%
    dplyr::filter(p.adjust < alpha) %>%
    dplyr::arrange(desc(NES))
  df
})


# Save
openxlsx::write.xlsx(
  markers_all,
  here("results/paper/tables/supplementary_table_myeloid_markers_all.xlsx")
)
openxlsx::write.xlsx(
  markers_slan_l,
  here("results/paper/tables/supplementary_table_myeloid_markers_slan.xlsx")
)
openxlsx::write.xlsx(
  markers_dc_l,
  here("results/paper/tables/supplementary_table_myeloid_markers_dc.xlsx")
)
openxlsx::write.xlsx(
  gsea_sorted,
  here("results/paper/tables/supplementary_table_myeloid_GO_slan.xlsx")
)
