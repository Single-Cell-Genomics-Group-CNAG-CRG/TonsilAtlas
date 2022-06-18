# This script finds all markers for all clusters of a specific cell type (one vs all)


# Load packages
library(tidyverse)
library(Seurat)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]


# Define paths
path_to_obj <- str_c(
  here::here("scRNA-seq/results/R_objects/level_3/"),
  cell_type,
  "/",
  cell_type,
  "_clustered_level_3.rds",
  sep = ""
)
path_to_tmp <- str_c(
  here::here("scRNA-seq/3-clustering/3-level_3/tmp/"),
  cell_type,
  sep = ""
)
dir.create(path_to_tmp, showWarnings = FALSE)
path_to_save_xlsx <- str_c(
  path_to_tmp,
  "/",
  cell_type,
  "_markers_level_3.xlsx"
)
path_to_save_rds <- str_c(
  path_to_tmp,
  "/",
  cell_type,
  "_markers_level_3.rds"
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
markers_dfs <- markers_dfs[as.character(sort(as.numeric(names(markers_dfs))))]
openxlsx::write.xlsx(markers_dfs, path_to_save_xlsx)
