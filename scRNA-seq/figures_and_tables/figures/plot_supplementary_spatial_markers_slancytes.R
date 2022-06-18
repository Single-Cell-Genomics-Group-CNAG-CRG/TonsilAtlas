# Load packages
library(Seurat)
library(tidyverse)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
seurat <- readRDS(here(path_to_save_tonsil))


# Dot plot
goi <- c("MMP12", "C1QA", "SELENOP")
Idents(seurat) <- "annotation_20220215"
dot_plot <- DotPlot(seurat, features = goi) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1))


# Save
ggsave(
  filename =  here("results/paper/figures/supplementary_figure_myeloid_dot_plot_spatial_markers_slancytes.pdf"),
  plot = dot_plot,
  device = cairo_pdf,
  width = 15, 
  height = 29, 
  units = "cm"
)
