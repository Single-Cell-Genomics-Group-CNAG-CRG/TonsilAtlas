# This script plots the UMAP and markers (RNA and CITE-seq) for CD8, ILC and
# NK cells

# Load packages
library(Seurat)
library(Signac)
library(scCustomize)
library(tidyverse)
library(ggrastr)
library(pheatmap2)
library(patchwork)
library(ggpubr)
library(here)
library(glue)


# Source utils
source(here("scRNA-seq/bin/utils_figure3.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))
path_to_wd <- here()
path_save_heatmap_rna <- here("results/paper/figures/figure_3_cytotoxic_heatmap_rna.pdf")
path_to_scirpy_output <- here("data/raw_data_figures/scirpy_tcr_output.tsv")


# Read data
cd8 <- readRDS(path_to_save_cd8)
ilc_nk <- readRDS(path_to_save_ilc_nk)
seurat_cite <- readRDS(path_to_save_cite_cytotoxic)
seurat_cite <- subset(seurat_cite, subset = annotation_prob >= 0.6)
atac <- readRDS(path_to_save_atac_cytotoxic)
scirpy_df <- read_csv(path_to_scirpy_output, col_names = TRUE)


# UMAP
ilc_nk$UMAP_1_20220215 <- ilc_nk$UMAP_1_20220215 + 10
ilc_nk$UMAP_2_20220215 <- ilc_nk$UMAP_2_20220215 + 10
seurat_rna <- merge(x = cd8, y = ilc_nk)
seurat_rna$annotation_20220215[seurat_rna$annotation_20220215 == "CXCR6+ RM CD8 T"] <- "RM CD8 activated T"
seurat_rna$annotation_paper <- factor(
  seurat_rna$annotation_20220215,
  levels = names(colors_rna)
)
seurat_rna <- ScaleData(seurat_rna, features = rownames(seurat_rna))
Idents(seurat_rna) <- "annotation_paper"
umap_title <- format(ncol(seurat_rna), big.mark = ",", scientific = FALSE)
umap_annotation <- seurat_rna@meta.data %>%
  ggplot(aes(UMAP_1_20220215, UMAP_2_20220215, color = annotation_paper)) +
  geom_point(shape = ".", alpha = 0.9) +
  scale_color_manual(values = colors_rna, breaks = names(colors_rna)) +
  labs(title = str_c(umap_title, "cells", sep = " "), x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  coord_fixed() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
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
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2)))
umap_annotation <- rasterize(umap_annotation, dpi = 300)


# Violin plot key cytotoxic markers
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
cell_types <- levels(seurat_rna$annotation_paper)
mycolors <- list(cell_type = c(colors_rna))
annotation_col <- data.frame(cell_type = cell_types) 
rownames(annotation_col) <- cell_types


# Scale per row between 0 and 1 
input_mat <- apply(avgexpr_mat, 1, function(x) (x - min(x)) / diff(range(x)))


# Plot heatmap
in2mm <- 25.4
pdf(path_save_heatmap_rna, width = (141 / in2mm), height = (80 / in2mm), paper = "special")
pheatmap2(
  input_mat,
  annotation_row = annotation_col,
  annotation_colors = mycolors,
  annotation_names_row = FALSE,
  annotation_legend = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE, 
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = c(8, 12, 15, 18),
  gaps_col = c(12, 29, 31, 35, 39),
  fontsize_row = 6
)
# Add black border to color boxes
# grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
# grid.gedit("col_annotation", gp = gpar(col = "black"))
dev.off()


# CITE-seq
DefaultAssay(seurat_cite) <- "ADT"
seurat_cite@reductions$wnn.umap@cell.embeddings[, "wnnUMAP_1"] <- seurat_cite$UMAP_1_RNA_based
seurat_cite@reductions$wnn.umap@cell.embeddings[, "wnnUMAP_2"] <- seurat_cite$UMAP_2_RNA_based
goi_cite <- c("CD45RA", "CD45RO", "CD99", "CD28", "CD38", "CD279-(PD-1)",
              "CD278-(ICOS)", "CD336-(NKp44)")
goi_cite2 <- c("CD103-(Integrin-alphaE)", "CD54", "CD161",
               "CD56-(NCAM)")
cite_umaps <- map(goi_cite, function(x) {
  p <- FeaturePlot(
    seurat_cite,
    x,
    reduction = "wnn.umap",
    pt.size = 0.01,
    order = FALSE
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
  rasterize(p, dpi = 300)
})
legend_umap_cite <- as_ggplot(get_legend(cite_umaps[[1]]))
cite_umaps <- purrr::map(cite_umaps, \(p) {
  p + theme(legend.position = "none")
})
cite_umaps2 <- map(goi_cite2, function(x) {
  p <- FeaturePlot(
    seurat_cite,
    x,
    reduction = "wnn.umap",
    pt.size = 0.01,
    order = TRUE
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
  rasterize(p, dpi = 300)
})
cite_umaps2 <- purrr::map(cite_umaps2, \(p) {
  p + theme(legend.position = "none")
})
cite_umaps <- c(cite_umaps, cite_umaps2)
cite_umaps <- wrap_plots(cite_umaps, ncol = 3)


# scTCR-seq
sel_cols <- c("barcode", "clonotype", "clonotype_size", "clonal_expansion_flag",
              "clonal_expansion")
scirpy_df <- scirpy_df[, sel_cols]
seurat_cite@meta.data <- left_join(
  seurat_cite@meta.data,
  scirpy_df,
  by = "barcode"
)
is_expanded <- 
  seurat_cite$clonal_expansion == ">= 3" &
  !is.na(seurat_cite$clonal_expansion)
selected_cells <- colnames(seurat_cite)[is_expanded]
colors_tcr <- c("#c3c3c3", "#f0b635")
umap_tcr <- DimPlot(
  seurat_cite,
  reduction = "wnn.umap",
  cells.highlight = selected_cells,
  pt.size = 0.1,
  sizes.highlight = 1
)
umap_tcr <- umap_tcr +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  ) +
  theme(
    legend.position = c(0, 0.75),
    legend.text = element_text(size = 6),
    plot.title = element_blank(),
    panel.border = element_blank(),
  ) +
  coord_fixed() +
  scale_color_manual(
    values = colors_tcr,
    labels = c("<3 cells/clonotype", ">=3 cells/clonotype")
  )
umap_tcr <- rasterize(umap_tcr, dpi = 300)
barplot_tcr_df <- seurat_cite@meta.data %>%
  mutate(annotation = factor(
    annotation_20220215,
    levels = names(colors_rna),
  )) %>% 
  dplyr::filter(clonal_expansion == ">= 3") %>%
  group_by(clonotype, annotation, .drop = FALSE) %>%
  summarise(n_cells = n())
levels_clonotype <- barplot_tcr_df %>%
  group_by(clonotype) %>%
  summarise(total_cells = sum(n_cells)) %>%
  arrange(desc(total_cells)) %>%
  pull(clonotype)
barplot_tcr_df$clonotype <- factor(barplot_tcr_df$clonotype, levels_clonotype)
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
fig_tcr <- umap_tcr | barplot_tcr_gg


# ATAC
DefaultAssay(atac) <- "chromvar"
moi_atac <- c("MA0652.1", "MA0624.1") # motif of interest
names(moi_atac) <- c("IRF8", "NFATC1")
atac@reductions$umap@cell.embeddings[, "UMAP_1"] <- atac$UMAP_1_RNA_based
atac@reductions$umap@cell.embeddings[, "UMAP_2"] <- atac$UMAP_2_RNA_based
umaps_accessibility <- purrr::map(moi_atac, function(x) {
  p <- FeaturePlot(
    atac,
    min.cutoff = "q5",
    max.cutoff = "q95",
    raster = FALSE,
    features = x,
    pt.size = 0.01,
    order = TRUE
  ) +
    scale_colour_gradient(low = "lightgray", high = "green4") +
    labs(x = "UMAP1", y = "UMAP2", color = "Accessibility") +
    theme_classic() +
    theme(
      legend.title = element_text(size = 6),
      legend.position = "bottom",
      legend.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  p
})
legend_umap_acc <- as_ggplot(get_legend(umaps_accessibility[[1]]))
umaps_accessibility[[1]] <- umaps_accessibility[[1]] +
  ggtitle("IRF8") +
  theme(plot.title = element_text(size = 6, hjust = 0.5))
umaps_accessibility[[2]] <- umaps_accessibility[[2]] +
  ggtitle("NFATC1") +
  theme(plot.title = element_text(size = 6, hjust = 0.5))
umaps_accessibility[[1]] <- rasterize(umaps_accessibility[[1]], dpi = 300)
umaps_accessibility[[2]] <- rasterize(umaps_accessibility[[2]], dpi = 300)
umaps2 <- (umaps_accessibility[[1]] / umaps_accessibility[[2]]) & NoLegend()
motif_plots <- map(moi_atac, function(x) {
  p <- MotifPlot(atac, x, assay = "peaks_redefined")
  p <- p +
    theme(
      plot.title = element_text(size = 5),
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 5)
    )
  p
})
motif_plots_arr <- wrap_plots(motif_plots, ncol = 1)

# Save
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_umap.pdf"),
  plot = umap_annotation,
  height = 8.5,
  width = 9,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_vln_markers.pdf"),
  plot = vlns,
  height = 15,
  width = 11.5,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_umaps_cite.pdf"),
  plot = cite_umaps,
  height = 16,
  width = 10,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_umaps_cite_legend.pdf"),
  plot = legend_umap_cite,
  height = 4,
  width = 7,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_atac_motifs.pdf"),
  plot = umaps2,
  height = 9,
  width = 6,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_atac_motifs_legend.pdf"),
  plot = legend_umap_acc,
  height = 4,
  width = 7,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_atac_motif_logoss.pdf"),
  plot = motif_plots_arr,
  height = 9,
  width = 6,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_tcr_umap.pdf"),
  plot = umap_tcr,
  height = 9,
  width = 6,
  units = "cm"
)
ggsave(
  filename = glue("{path_to_wd}/results/paper/figures/figure_3_cytotoxic_tcr_barplot.pdf"),
  plot = barplot_tcr_gg,
  height = 5.5,
  width = 5.5,
  units = "cm"
)
