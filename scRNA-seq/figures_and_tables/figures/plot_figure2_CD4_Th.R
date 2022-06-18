# This script plots the UMAP and markers for the different subsets of Th (main
# figure 2)


# Load packages
library(Seurat)
library(Nebulosa)
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(here)


# Source utils
source(here("scRNA-seq/bin/utils_figure2.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
seurat_th <- readRDS(path_to_save_th)


# Plot Th expression and markers
umap_th <- DimPlot(seurat_th, pt.size = 0.1, cols = colors_th)
umap_th <- umap_th +
  theme_nothing2() +
  theme(
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    legend.spacing.x = unit(0.01, "cm")
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 1.5)))
umaps_th_markers <- plot_density(
  seurat_th,
  goi_th_main,
  reduction = "umap",
  size = 0.1,
  combine = FALSE
)
umaps_th_markers <- map(umaps_th_markers, function(p) {
  p +   scale_color_viridis_c(option = "magma") &
    NoLegend() &
    theme(
      plot.title = element_text(size = 5, hjust = 0.5),
      axis.text = element_text(size = 0),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
})
umaps_th_markers <- ggarrange(plotlist = umaps_th_markers, ncol = 4, nrow = 2)
umaps_th_arr <- umap_th / umaps_th_markers


# Supplementary
# umaps_th_supp <- plot_density(
#   seurat_th,
#   goi_th_supp,
#   reduction = "umap",
#   size = 0.1,
#   combine = FALSE
# )
# umaps_th_supp <- map(umaps_th_supp, function(p) {
#   p + scale_color_viridis_c(option = "magma") &
#     NoLegend() &
#     theme(
#       plot.title = element_text(size = 5, hjust = 0.5),
#       axis.text = element_text(size = 0),
#       axis.line = element_blank(),
#       axis.title = element_blank(),
#       axis.ticks = element_blank()
#     )
# })
# umaps_th_supp <- ggarrange(plotlist = umaps_th_supp, ncol = 5, nrow = 2)


# Save
ggsave(
  filename = path_to_save_th_main,
  plot = umaps_th_arr,
  width = 7,
  height = 9.5,
  units = "cm"
)
ggsave(
  filename = path_to_save_th_supp,
  plot = umaps_th_supp,
  width = 21,
  height = 9.5,
  units = "cm"
)
