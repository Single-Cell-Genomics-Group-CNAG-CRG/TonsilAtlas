# This script finds all markers for all clusters of a specific cell type (one vs all)


# Load packages
library(tidyverse)
library(Seurat)


# Define paths
path_to_obj <- here::here("scRNA-seq/results/R_objects/level_4/CD4_T/CD4_T_integrated_level_4.rds") 
path_to_tmp <- here::here("scRNA-seq/3-clustering/4-level_4/CD4_T/tmp/")
dir.create(path_to_tmp, showWarnings = FALSE, recursive = TRUE)
path_to_save_xlsx <- str_c(
  path_to_tmp,
  "/CD4_T_markers_level_4.xlsx"
)
path_to_save_rds <- str_c(
  path_to_tmp,
  "/CD4_T_markers_level_4.rds"
)


# Load data
seurat <- readRDS(path_to_obj)


# Find markers
markers <- FindAllMarkers(
  seurat,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  only.pos = TRUE,
  verbose = TRUE
)
markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_logFC), group_by = TRUE) %>%
  ungroup()


# Save
saveRDS(markers, path_to_save_rds)
markers_dfs <- purrr::map(unique(markers$cluster), function(x) {
  df <- markers[markers$cluster == x, ]
  df <- df[, c(7, 1:6)]
  df
})
names(markers_dfs) <- unique(markers$cluster)
openxlsx::write.xlsx(markers_dfs, path_to_save_xlsx)
