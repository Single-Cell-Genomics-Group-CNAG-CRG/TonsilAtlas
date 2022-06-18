# This script plots the UMAP and markers (RNA and CITE-seq) for main figure 2
# (CD4 T cells)


# Load packages
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(grid)
library(patchwork)
library(pheatmap2)
library(here)


# Source utils
source(here("scRNA-seq/bin/utils_figure2.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
seurat_rna <- readRDS(path_to_save_cd4)
seurat_cite <- readRDS(path_to_save_cite_cd4)


# UMAP
seurat_rna$annotation_20220215 <- unfactor(seurat_rna$annotation_20220215)
seurat_rna$annotation_20220215[seurat_rna$annotation_20220215 == "non-GC-Tf-regs"] <- "Eff-Tregs-IL32"
seurat_rna$annotation_20220215[seurat_rna$annotation_20220215 == "GC-Tf-regs"] <- "Tfr"
seurat_rna$annotation_paper <- factor(
  seurat_rna$annotation_20220215,
  levels = names(colors_rna)
)
Idents(seurat_rna) <- "annotation_paper"
umap_title <- format(ncol(seurat_rna), big.mark = ",", scientific = FALSE)
umap_annotation <- Embeddings(seurat_rna, "umap") %>%
  as.data.frame() %>%
  mutate(annotation_paper = seurat_rna$annotation_paper) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = annotation_paper)) +
    geom_point(shape = ".", alpha = 0.9) +
    scale_color_manual(values = colors_rna, breaks = names(colors_rna)) +
    labs(title = str_c(umap_title, "cells", sep = " "), x = "UMAP1", y = "UMAP2") +
    theme_classic() +
    coord_fixed() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.spacing = unit(0, "cm"),
      legend.box.spacing = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.key.height = unit(0, "cm"),
      legend.key.width = unit(0, "cm"),
      plot.title = element_text(hjust = 0.5, size = 7),
      axis.title = element_text(size = 6),
      legend.text = element_text(size = 5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    ) +
    guides(
      color = guide_legend(
        override.aes = list(shape = 16, size = 2),
        nrow = 7,
        byrow = FALSE
      )
    )
umap_annotation <- rasterize(umap_annotation, dpi = 300)


# HEATMAP RNA
# Compute average expression per cluster
avgexpr_mat <- AverageExpression(
  seurat_rna,
  features = goi_rna,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "annotation_paper",
  slot = "data"
)$RNA


# Define annotation column
cell_types <- levels(seurat_rna$annotation_20220215)
mycolors <- list(cell_type = c(colors_rna))
annotation_col <- data.frame(cell_type = cell_types) 
rownames(annotation_col) <- cell_types


# Scale per row between 0 and 1 
input_mat <- t(apply(avgexpr_mat, 1, function(x) (x - min(x)) / diff(range(x))))


# Plot heatmap
in2mm <- 25.4
pdf(path_save_heatmap_rna, width = (70 / in2mm), height = (125 / in2mm), paper = "special")
pheatmap2(
  input_mat,
  annotation_col = annotation_col,
  annotation_colors = mycolors,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE, 
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = c(1, 3, 8),
  fontsize_row = 6
)
# Add black border to color boxes
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
grid.gedit("col_annotation", gp = gpar(col = "black"))
dev.off()


# Markers CITE-seq
DefaultAssay(seurat_cite) <- "ADT"
seurat_cite@reductions$wnn.umap@cell.embeddings[, "wnnUMAP_1"] <- seurat_cite$UMAP_1_20220215
seurat_cite@reductions$wnn.umap@cell.embeddings[, "wnnUMAP_2"] <- seurat_cite$UMAP_2_20220215
cite_umaps <- map(goi_cite, function(x) {
  p <- FeaturePlot(
    seurat_cite,
    x,
    reduction = "wnn.umap",
    pt.size = 0.001,
    min.cutoff = "q1",
    max.cutoff = "q99",
    order = FALSE,
    slot = "data"
  ) +
    scale_color_viridis_c(option = "magma") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 6, face = "plain", hjust = 0.5),
      legend.title = element_text(size = 6, face = "plain"),
      legend.position = "bottom",
      legend.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  p$layers[[1]]$aes_params$shape <- "."
  p$layers[[1]]$aes_params$size <- NULL
  rasterize(p, dpi = 720)
})
legend_umap_cite <- as_ggplot(get_legend(cite_umaps[[1]]))
cite_umaps <- purrr::map(cite_umaps, \(p) {
  p + theme(legend.position = "none")
})
cite_umaps_arr <- wrap_plots(cite_umaps, ncol = 3, nrow = 3)


# UMAP groups
colors_groups <- c(
  "Naive" = "#dedede",
  "CM" = "#7d7d7d",
  "Non-Tfh" = "#fa7f70",
  "Tfh" = "#3f68e1"
)
umap_new_groups <- Embeddings(seurat_rna, "umap") %>%
  as.data.frame() %>%
  mutate(annotation_paper = seurat_rna$annotation_paper) %>%
  mutate(new_groups = case_when(
    annotation_paper == "Naive" ~ "Naive",
    annotation_paper %in% c("CM PreTfh", "CM Pre-non-Tfh") ~ "CM",
    annotation_paper %in% c("T-Trans-Mem", "T-helper", "Eff-Tregs", "GC-Tf-regs", "non-GC-Tf-regs") ~ "Non-Tfh",
    annotation_paper %in% c("Tfh-LZ-GC", "GC-Tfh-SAP", "GC-Tfh-0X40", "Tfh-Mem", "Tfh T:B border") ~ "Tfh",
  )) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = new_groups)) +
  geom_point(shape = ".", alpha = 0.9) +
  scale_color_manual(values = colors_groups, breaks = names(colors_groups)) +
  theme_classic() +
  coord_fixed() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
umap_new_groups <- rasterize(umap_new_groups, dpi = 300)


# Save
ggsave(
  plot = umap_annotation,
  filename = path_save_umap_fig2,
  width = 8.5,
  height = 8.5,
  units = "cm"
)
# ggsave(
#   plot = heatmap_rna,
#   filename = path_save_heatmap_rna,
#   width = 7.5,
#   height = 14,
#   units = "cm"
# )
# ggsave(
#   plot = legend_heatmap,
#   filename = path_save_heatmap_rna_leg,
#   width = 15,
#   height = 3,
#   units = "cm"
# )
ggsave(
  plot = cite_umaps_arr,
  filename = path_save_umaps_cite,
  width = 7,
  height = 7,
  units = "cm"
)
ggsave(
  plot = legend_dot_plot,
  filename = path_save_dotplot_cite_leg,
  width = 15,
  height = 3,
  units = "cm"
)
ggsave(
  plot = umap_new_groups,
  filename = here("results/paper/figures/figure_2_umap_main_groups_rna.pdf"),
  width = 4,
  height = 4,
  units = "cm"
)
