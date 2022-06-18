# This script creates the pdfs with the relevant plots to share with the
# annotation team


# Load packages
library(png)
library(grid)
library(gridExtra)
library(tidyverse)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]


# Define paths
path_to_tmp <- str_c(
  here::here("scRNA-seq/3-clustering/3-level_3/tmp/"),
  cell_type,
  "/",
  sep = ""
)
path_to_umap_pre_freeze <- str_c(path_to_tmp, cell_type, "_umap_annotation_pre_freeze.png")
path_to_barplot_pre_freeze <- str_c(path_to_tmp, cell_type, "_barplot_annotation_pre_freeze.png")
path_to_save_pre_freeze <- str_c(path_to_tmp, cell_type, "_annotation_pre_freeze.pdf")

path_to_umap_king <- str_c(path_to_tmp, cell_type, "_umap_annotation_king_et_al_level_3.png")
path_to_barplot_king <- str_c(path_to_tmp, cell_type, "_barplot_annotation_king_et_al_level_3.png")
path_to_save_annotation_king <- str_c(path_to_tmp, cell_type, "_annotation_king.pdf")

path_to_umap_level_1 <- str_c(path_to_tmp, cell_type, "_umap_clusters_level_1.png")
path_to_umap_level_3 <- str_c(path_to_tmp, cell_type, "_umap_clusters_level_3.png")
path_to_save_umaps_clusters <- str_c(path_to_tmp, cell_type, "_umaps_clusters.pdf")



# Read PNG
umap_pre_freeze <- rasterGrob(readPNG(path_to_umap_pre_freeze, native = FALSE), interpolate = FALSE)
barplot_pre_freeze <- rasterGrob(readPNG(path_to_barplot_pre_freeze, native = FALSE), interpolate = FALSE)

umap_king <- rasterGrob(readPNG(path_to_umap_king, native = FALSE), interpolate = FALSE)
barplot_king <- rasterGrob(readPNG(path_to_barplot_king, native = FALSE), interpolate = FALSE)

umap_level_1 <- rasterGrob(readPNG(path_to_umap_level_1, native = FALSE), interpolate = FALSE)
umap_level_3 <- rasterGrob(readPNG(path_to_umap_level_3, native = FALSE), interpolate = FALSE)


# Create and save pdfs
pdf(path_to_save_pre_freeze)
do.call(grid.arrange, list(umap_pre_freeze))
do.call(grid.arrange, list(barplot_pre_freeze))
dev.off()

pdf(path_to_save_annotation_king)
do.call(grid.arrange, list(umap_king))
do.call(grid.arrange, list(barplot_king))
dev.off()

pdf(path_to_save_umaps_clusters)
do.call(grid.arrange, list(umap_level_3))
do.call(grid.arrange, list(umap_level_1))
dev.off()
