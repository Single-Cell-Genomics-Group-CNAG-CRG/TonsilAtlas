# This script plots the results of scTCR-seq analysis in figure 2


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(here)
library(ggrastr)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure2.R"))


# Read data
seurat_cite <- readRDS(path_to_save_cite_cd4)


# Plot
seurat_cite@reductions$wnn.umap@cell.embeddings[, 1] <- seurat_cite$UMAP_1_20220215
seurat_cite@reductions$wnn.umap@cell.embeddings[, 2] <- seurat_cite$UMAP_2_20220215
is_expanded <- seurat_cite$clonal_expansion == ">= 3"
selected_cells <- colnames(seurat_cite)[is_expanded]
umap_tcr <- DimPlot(
  seurat_cite,
  reduction = "wnn.umap",
  cells.highlight = selected_cells,
  pt.size = 0.1,
  sizes.highlight = 1
)
umap_tcr <- rasterize(umap_tcr, dpi = 300)
umap_tcr <- umap_tcr +
  theme_nothing2() +
  theme(
    legend.position = c(0.65, 0.25),
    legend.text = element_text(size = 6),
    plot.title = element_blank(),
    panel.border = element_blank(),
  ) +
  coord_fixed() +
  scale_color_manual(
    values = colors_tcr,
    labels = c("<3 cells/clonotype", ">=3 cells/clonotype")
  )
barplot_tcr_df <- seurat_cite@meta.data %>%
  mutate(annotation = factor(
    annotation_figure_2,
    levels = names(colors_rna),
  )) %>% 
  dplyr::filter(clonal_expansion == ">= 3") %>%
  group_by(clonotype, annotation, .drop = FALSE) %>%
  summarise(n_cells = n())
barplot_tcr_df$clonotype <- factor(
  barplot_tcr_df$clonotype,
  levels = c("2493_TCR", "2072_TCR", "3573_TCR", "3699_TCR", "862_TCR",
             "879_TCR", "901_TCR")
)
barplot_tcr_gg <- barplot_tcr_df %>%
  ggplot(aes(clonotype, n_cells, fill = annotation)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = colors_rna) +
  scale_x_discrete(
    breaks = levels(barplot_tcr_df$clonotype),
    labels = 1:length(levels(barplot_tcr_df$clonotype))
  ) +
  NoLegend() +
  ylab("Number of cells") +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7)
  )


# Arrange
tcr <- umap_tcr / barplot_tcr_gg


# Save
ggsave(
  filename = path_to_save_figure_2_scTCR_seq,
  tcr,
  width = 3.5,
  height = 7,
  units = "cm"
)
