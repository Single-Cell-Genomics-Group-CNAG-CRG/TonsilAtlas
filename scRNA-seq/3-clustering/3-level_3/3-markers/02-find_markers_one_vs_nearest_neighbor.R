# This script finds all the markers between one cluster and its most similar
# cluster
# https://dplyr.tidyverse.org/reference/across.html


# Load packages
library(Seurat)
library(tidyverse)


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
path_to_save_xlsx <- str_c(
  path_to_tmp,
  "/",
  cell_type,
  "_markers_level_3_one_vs_nearest_neighbor.xlsx"
)
path_to_save_rds <- str_c(
  path_to_tmp,
  "/",
  cell_type,
  "_markers_level_3_one_vs_nearest_neighbor.rds"
)
path_to_annotation <- str_c(here::here("annotation/level_3/"), cell_type)
path_to_save_xlsx_annotation <- str_c(
  path_to_annotation,
  "/",
  cell_type,
  "_markers_level_3_one_vs_nearest_neighbor.xlsx"
)


# Load data
seurat <- readRDS(path_to_obj)


# Find nearest neighbor for each cluster
pseudo_bulk_harmony_df <- seurat@reductions$harmony@cell.embeddings[, 1:30] %>%
  as.data.frame() %>%
  mutate(seurat_clusters = seurat$seurat_clusters) %>%
  group_by(seurat_clusters) %>%
  summarise(across(starts_with("harmony_"), ~mean(.x, na.rm = TRUE)))
pseudo_bulk_harmony_mat <- pseudo_bulk_harmony_df %>%
  select(starts_with("harmony_")) %>%
  as.matrix()
rownames(pseudo_bulk_harmony_mat) <- pseudo_bulk_harmony_df$seurat_clusters
dist_mat <- as.matrix(dist(pseudo_bulk_harmony_mat, upper = TRUE, method = "euclidean"))
diag(dist_mat) <- NA
n_cells_by_cluster <- table(seurat$seurat_clusters)
min_cells <- ifelse(cell_type != "epithelial", 150, 25)
nearest_neighbors <- purrr::map(rownames(dist_mat), function(x) {
  dists <- dist_mat[x, ]
  nn <- names(which.min(dists))
  if (n_cells_by_cluster[nn] > min_cells) {
    return(nn)
  } else {
    nn2 <- names(which.min(dists[-which(nn == names(dists))]))
    return(c(nn, nn2))
  }
})
names(nearest_neighbors) <- rownames(dist_mat)


# Find Markers
markers_dfs <- purrr::map(names(nearest_neighbors), function(x) {
  df <- FindMarkers(
    seurat,
    ident.1 = x,
    ident.2 = nearest_neighbors[[x]],
    logfc.threshold = 0.5,
    test.use = "wilcox",
    only.pos = FALSE,
    verbose = TRUE
  )
  df <- df %>%
    rownames_to_column(var = "gene") %>% 
    arrange(desc(avg_logFC))
  df
})
names_dfs <- purrr::map2_chr(names(nearest_neighbors), nearest_neighbors, function(cluster, nns) {
  nns <- str_c(nns, collapse = ";")
  out <- str_c(cluster, nns, sep = "vs")
  out
})
names(markers_dfs) <- names_dfs


# Save
saveRDS(markers_dfs, path_to_save_rds)
openxlsx::write.xlsx(markers_dfs, path_to_save_xlsx)
openxlsx::write.xlsx(markers_dfs, path_to_save_xlsx_annotation)


