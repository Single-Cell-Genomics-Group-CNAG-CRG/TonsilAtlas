# This script plots all the results concerning PRDM1 and BCL6 from figure 2:
# - Regulons
# - scATAC-seq


# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(ggrastr)
library(ggplotify)
library(ggpubr)
library(pheatmap2)
library(patchwork)
library(GenomicRanges)
library(here)


# Source utils
source(here("scRNA-seq/bin/utils_figure2.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
seurat_rna <- readRDS(path_to_save_cd4)
seurat_atac <- readRDS(path_to_save_atac_cd4)
auc_mtx <-  read.csv(path_to_auc_mtx, row.names = 1, check.names = FALSE)
rna_df <- read_delim(path_to_mat_rna_bcl6, delim = " ", skip = 1, col_names = FALSE)
atac_df <- read_delim(path_to_mat_atac_bcl6, delim = " ", skip = 1, col_names = FALSE)


# Plot activity PRDM1 and BCL6
filt <- colnames(seurat_rna)[colnames(seurat_rna) %in% rownames(auc_mtx)]
auc_mtx <- auc_mtx[filt, c("BCL6(+)", "PRDM1(+)")]
seurat_rna$BCL6_activity <- NA
seurat_rna$PRDM1_activity <- NA
seurat_rna$BCL6_activity[rownames(auc_mtx)] <- auc_mtx[, "BCL6(+)"]
seurat_rna$PRDM1_activity[rownames(auc_mtx)] <- auc_mtx[, "PRDM1(+)"]
vars1 <- c("BCL6_activity", "PRDM1_activity")
seurat_rna_sub <- subset(
  seurat_rna,
  cells = colnames(seurat_rna)[!is.na(seurat_rna$BCL6_activity)]
)
umaps_activities <- map(vars1, function(x) {
  p <- FeaturePlot(
    seurat_rna_sub,
    raster = FALSE,
    features = x,
    pt.size = 0.01,
    order = TRUE
  ) +
    scale_colour_gradient(low = "lightgrey", high = "red") +
    labs(x = "UMAP1", y = "UMAP2", color = "TF activity") +
    theme_classic() +
    coord_fixed() +
    theme(
      legend.title = element_text(size = 6),
      legend.position = "bottom",
      legend.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
   rasterize(p, dpi = 300)
})
legend_umap_act <- as_ggplot(get_legend(umaps_activities[[1]]))
umaps_activities[[1]] <- umaps_activities[[1]] +
  ggtitle("BCL6") +
  theme(legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5))
umaps_activities[[2]] <- umaps_activities[[2]] +
  ggtitle("PRDM1") +
  theme(legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5))
umaps <- umaps_activities[[1]] / umaps_activities[[2]]


# Expression vs accessibility
atac_mat <- as.matrix(atac_df[, 2:ncol(atac_df)])
rownames(atac_mat) <- atac_df$X1
colnames(atac_mat) <- c("CM", "Naive", "Non-Tfh", "Tfh")
new_order <- c("Naive", "CM", "Tfh", "Non-Tfh")
atac_mat <- atac_mat[, new_order]
rna_mat <- as.matrix(rna_df[, 2:ncol(rna_df)])
rownames(rna_mat) <- rna_df$X1
colnames(rna_mat) <- c("CM", "Naive", "Non-Tfh", "Tfh")
rna_mat <- rna_mat[, new_order]
in2mm <- 25.4
pdf(path_save_heatmaps_tfh_rna, width = (65 / in2mm), height = (16 / in2mm), paper = "special")
pheatmap2(
  t(rna_mat),
  scale = "column",
  annotation_names_row = FALSE,
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE
)
dev.off()
colfunc <- colorRampPalette(c("#4575B4", "#FFFFBF", "darkgreen"))
pdf(path_save_heatmaps_tfh_atac, width = (65 / in2mm), height = (16 / in2mm), paper = "special")
pheatmap2(
  t(atac_mat),
  scale = "column",
  annotation_names_row = FALSE,
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = FALSE,
  color = colfunc(100),
  show_rownames = FALSE,
  show_colnames = FALSE
)
dev.off()
pdf(path_save_heatmaps_tfh_atac_legend, width = (70 / in2mm), height = (18 / in2mm), paper = "special")
pheatmap2(
  t(atac_mat),
  scale = "column",
  annotation_names_row = FALSE,
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_col = 6,
  angle_col = "90",
  legend = TRUE,
  color = colfunc(100),
  show_rownames = FALSE,
  show_colnames = TRUE
)
dev.off()


# Coverage plot BCL6 superenhancer

# Restingir a 18760...
# Scale peaks size
# Canviar ordre: 1. genes, 2. coverage, 3. links (like genome browser). Eliminar peaks
# Inkscape: canviar size i posició gene labels
# Inkscape:  eliminar annotation label (Naive...), afegir "boles"
# Track: canviar a gris?
# region_highlight_enhancer <- GRanges(
#   seqnames = "chr3",
#   ranges = IRanges(start = 187650000, end = 188001000)
# )
tfh_levels <- c("Tfh T:B border", "Tfh-LZ-GC", "GC-Tfh-SAP", "GC-Tfh-0X40",
                "Tfh-Mem")
seurat_atac$annotation2 <- unfactor(seurat_atac$annotation_20220215)
seurat_atac$annotation2[!(seurat_atac$annotation2 %in% tfh_levels)] <- "Other CD4 T"
seurat_atac$annotation2 <- factor(
  seurat_atac$annotation2,
  levels = c("Other CD4 T", tfh_levels)
)
Idents(seurat_atac) <- "annotation2"
region_bcl6_links <- "chr3-187725000-188000857"
coverage_gg <- CoveragePlot(
  object = seurat_atac,
  # region.highlight = region_highlight_enhancer,
  group.by = "annotation2",
  region = region_bcl6_links,
  heights = c(16, 1, 1, 4)
) &
  scale_fill_manual(values = rep("black", length(colors_rna)), breaks = names(colors_rna))
# coverage_gg <- rasterize(coverage_gg, dpi = 300)
# coverage_gg[[1]] <- coverage_gg[[1]] +
#   scale_fill_manual(values = rep("gray37", length(colors_rna)), breaks = names(colors_rna))
# coverage_gg[[2]]$layers[[5]]$aes_params$size <- 2
# coverage_gg[[4]] <- coverage_gg[[4]] +
#   theme(
#     axis.title.x.bottom = element_text(size = 7),
#     axis.text.x = element_text(size = 6),
#     legend.title = element_text(size = 6.5),
#     legend.text = element_text(size = 6)
#   )
# coverage_gg

# Plot accessibility PRDM1 and BCL6
vars2 <- c("BCL6_gene1", "BCL6_enhancer1")
seurat_atac@reductions$umap@cell.embeddings[, "UMAP_1"] <- seurat_atac$UMAP_1_RNA_based
seurat_atac@reductions$umap@cell.embeddings[, "UMAP_2"] <- seurat_atac$UMAP_2_RNA_based
umaps_accessibility <- purrr::map(vars2, function(x) {
  p <- FeaturePlot(
    seurat_atac,
    min.cutoff = "q5",
    max.cutoff = "q95",
    raster = FALSE,
    features = x,
    pt.size = 0.01,
    order = TRUE
  ) +
    scale_colour_gradient(low = "lightgray", high = "green4") +
    labs(x = "UMAP1", y = "UMAP2", color = "Accessibility Score") +
      theme_classic() +
      coord_fixed() +
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
  ggtitle("BCL6 ATAC Score (promoter)") +
  theme(legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5))
umaps_accessibility[[2]] <- umaps_accessibility[[2]] +
  ggtitle("BCL6 ATAC Score (enhancer)") +
  theme(legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5))
umaps_accessibility[[1]] <- rasterize(umaps_accessibility[[1]], dpi = 300)
umaps_accessibility[[2]] <- rasterize(umaps_accessibility[[2]], dpi = 300)
umaps2 <- umaps_accessibility[[1]] / umaps_accessibility[[2]]
umaps2
# Add number of cells as annotation
# Change green to same as heatmap


# Save
ggsave(
  plot = umaps,
  filename = path_save_umaps_act,
  width = 3.5,
  height = 7,
  units = "cm"
)
ggsave(
  plot = legend_umap_act,
  filename = path_save_umaps_act_leg,
  width = 3.5,
  height = 7,
  units = "cm"
)
# ggsave(
#   plot = heatmap_tfh,
#   filename = path_save_heatmaps_tfh,
#   width = 6,
#   height = 7,
#   units = "cm"
# )
ggsave(
  plot = umaps2,
  filename = path_save_umaps_acc,
  width = 3.5,
  height = 7,
  units = "cm"
)
ggsave(
  plot = legend_umap_acc,
  filename = path_save_umaps_acc_leg,
  width = 3.5,
  height = 7,
  units = "cm"
)

# Save components coverage plot
# coverage_gg_main <- rasterize(coverage_gg[[1]], dpi = 300)
# coverage_gg_genes <- rasterize(coverage_gg[[2]], dpi = 300)

# test <- (coverage_gg[[1]] / coverage_gg[[2]] / coverage_gg[[4]]) +
# plot_layout(heights = c(1, 0.2, 0.2))
ggsave(
  plot = coverage_gg[[1]][[1]],
  filename = here("results/paper/figures/figure_2_coverage_plot_tfh_main.pdf"),
  width = 9.5,
  height = 6,
  units = "cm"
)
ggsave(
  plot = coverage_gg[[1]][[2]],
  filename = here("results/paper/figures/figure_2_coverage_plot_tfh_genes.pdf"),
  width = 9.5,
  height = 3,
  units = "cm"
)
ggsave(
  plot = coverage_gg[[1]][[3]],
  filename = here("results/paper/figures/figure_2_coverage_plot_tfh_links.pdf"),
  width = 9.5,
  height = 3,
  units = "cm"
)
