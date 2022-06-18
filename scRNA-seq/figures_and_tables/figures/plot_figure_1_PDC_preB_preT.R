# This script plots the supplementary figure containing:
# 1. UMAP with PDC, precursor T cells and precursor B cells
# 2. Dot plot with key markers


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(here)


# Define parameters
path_to_obj <- here("data/raw_data_figures/seurat_obj_pdc_preB_preT_fig1.rds")
cols <- c(
  "PDC" = "#f032e6",
  "IFN1+ PDC" = "#510791",
  "preB" = "#4363d8",
  "preT" = "#43d8c9"
)


# Read data
seurat <- readRDS(path_to_obj)


# UMAP
# seurat$annotation_20220215[seurat$annotation_20220215 == "IFN1+ PDC"] <- "PDC"
(umap_gg <- Embeddings(seurat, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = seurat$annotation_20220215) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = annotation)) +
    geom_point(size = 0.5, alpha = 0.75) +
    scale_color_manual(breaks = names(cols), values = cols) +
    labs(x = "UMAP1", y = "UMAP2") +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = c(0.85, 0.25),
      legend.text = element_text(size = 5),
      axis.text = element_blank(),
      axis.title = element_text(size = 5),
      axis.ticks = element_blank()
    )) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# Dot Plot
seurat$annotation_20220215 <- factor(
  seurat$annotation_20220215,
  levels = rev(c("PDC", "IFN1+ PDC", "preB", "preT"))
)
Idents(seurat) <- "annotation_20220215"
features <-  rev(c("IRF8", "IL3RA", "LILRA4", "IFI44", "IFI6", "CD79B", "CD19",
                   "CD3G", "CD8A", "PAX5", "RAG1", "RAG2", "DNTT",
                   "CD1A", "CD1B", "CD1C", "CD1E"))
dot_plot <- DotPlot(seurat, features = rev(features)) +
  # coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  scale_size(range = c(0, 4)) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 1),
    legend.position = "bottom")
legend_dot_plot <- as_ggplot(get_legend(dot_plot))
dot_plot <- dot_plot + theme(legend.position = "none")


# Save
ggsave(
  filename = here("results/paper/figures/figure_1_PDC_preB_preT_umap.pdf"),
  plot = umap_gg,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_1_PDC_preB_preT_dotplot.pdf"),
  plot = dot_plot,
  width = 6.5,
  height = 4.5,
  units = "cm"
)

ggsave(
  filename = here("results/paper/figures/figure_1_PDC_preB_preT_dotplot_legend.pdf"),
  plot = legend_dot_plot,
  width = 14,
  height = 4.5,
  units = "cm"
)
