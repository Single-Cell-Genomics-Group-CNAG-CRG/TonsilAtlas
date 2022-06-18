# This script plots figure 2


# Load packages
library(Seurat)
library(Signac)
library(GenomicRanges)
library(scCustomize)
library(ComplexHeatmap)
library(tidyverse)
library(ggpubr)
library(ggmap)
library(ggplotify)
library(ggrastr)
library(Nebulosa)
library(patchwork)
library(dittoSeq)
set.seed(123)


# Source utils
path_to_utils <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/bin/utils_figure2.R"
source(path_to_utils)


# Read data
seurat <- readRDS(path_to_rna)
seurat_cite <- readRDS(path_to_cite)
df_rna <- read_delim(path_to_mat_rna, delim = " ", skip = 1, col_names = FALSE)
df_atac <- read_delim(path_to_mat_atac, delim = " ", skip = 1, col_names = FALSE)
seurat_atac <- readRDS(path_to_atac)
seurat_th <- readRDS(path_to_th)
auc_mtx_all <-  read.csv(path_to_auc_mtx, row.names = 1, check.names = FALSE)
spata_obj <- readRDS(path_to_spata)
seurat_spatial <- readRDS(path_to_spatial)


# Plot UMAP
seurat$annotation_paper <- factor(
  seurat$annotation_paper,
  levels = names(colors_rna)
)
Idents(seurat) <- "annotation_paper"
umap_rna <- DimPlot(
  seurat,
  group.by = "annotation_paper",
  raster = FALSE,
  cols = colors_rna,
  pt.size = 0.005
) %>%
  rasterize(dpi = 300)
umap_rna <- umap_rna +
  theme_nothing2() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.spacing.x = unit(0.01, "cm"),
    plot.title = element_blank(),
    panel.border = element_blank()
  ) +
  coord_fixed() +
  guides(color = guide_legend(nrow = 6, ncol = 3, byrow = FALSE,
                              override.aes = list(size = 1.5)))
legend_umap <- as_ggplot(get_legend(umap_rna))
umap_rna <- umap_rna + theme(legend.position = "none")


# Plot heatmap with genes of interest
goi_rna <- unlist(goi_rna)
seurat_sub <- downsample_cells2(
  seurat,
  var = "annotation_paper",
  n_cells = 5000
)
seurat_sub <- ScaleData(seurat_sub, features = goi_rna)
heatmap_rna <- DoHeatmap(
  seurat_sub,
  features = goi_rna,
  group.by = "annotation_paper",
  slot = "scale.data",
  label = FALSE,
  raster = FALSE,
  group.colors = colors_rna,
) &
  theme(
    axis.text = element_text(size = 5),
    legend.position = "bottom",
    legend.title = element_text(size = 6),
    legend.text = element_blank,
    legend.spacing.x = unit(0.01, "cm")
  )
legend_heatmap <- as_ggplot(get_legend(heatmap_rna))
heatmap_rna <- heatmap_rna & no_legend()
heatmap_rna <- rasterize(heatmap_rna, dpi = 300)
empty_plot <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_blank() +
  theme_nothing()
heatmap_rna <- (heatmap_rna / empty_plot) +
  plot_layout(heights = c(0.975, 0.025))


# Plot CITE-seq
goi_cite <- unlist(goi_cite)
names(goi_cite) <- NULL
seurat_cite$annotation_paper <- factor(
  seurat_cite$CD4_annotation_level_5,
  levels = names(colors_rna)
)
dot_plot_cite <- DotPlot(
  seurat_cite,
  features = rev(goi_cite),
  group.by = "annotation_paper"
) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  scale_size_continuous(range = c(0.1, 4.5)) +
  coord_flip()
goi_cite_lables <- goi_cite
goi_cite_lables[goi_cite_lables == "CD127-(IL-7Ralpha)"] <- "CD127"
names(goi_cite_lables) <- goi_cite
dot_plot_cite <- dot_plot_cite +
  scale_x_discrete(labels = goi_cite_lables) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.spacing.x = unit(0.01, "cm"),
    legend.position = "bottom"
  )
legend_dot_plot <- as_ggplot(get_legend(dot_plot_cite))
dot_plot_cite <- dot_plot_cite & NoLegend()
tiles_df <- data.frame(annotation = names(colors_rna), x = 1:length(colors_rna), y = 1)
tiles_gg <-ggplot(tiles_df, aes(x, y, color = annotation)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors_rna) +
  theme_nothing() 
dot_plot_cite <- (dot_plot_cite / tiles_gg) +
  plot_layout(heights = c(0.975, 0.025))


# Plot activity PRDM1 and BCL6
auc_mtx_all <- auc_mtx_all[colnames(seurat), ]
auc_mtx_all_int_reg <- auc_mtx_all[, c("BCL6(+)", "PRDM1(+)")]
seurat$BCL6_activity <- NA
seurat$PRDM1_activity <- NA
seurat$BCL6_activity[rownames(auc_mtx_all_int_reg)] <- auc_mtx_all_int_reg[, "BCL6(+)"]
seurat$PRDM1_activity[rownames(auc_mtx_all_int_reg)] <- auc_mtx_all_int_reg[, "PRDM1(+)"]
seurat_rna <- subset(seurat, assay == "3P")
umaps_activities <- purrr::map(c("BCL6_activity", "PRDM1_activity"), function(x) {
  p <- FeaturePlot(
    seurat_rna,
    features = x,
    pt.size = 0.005,
    raster = FALSE
  ) +
    scale_color_viridis_c(option = "magma") +
    coord_fixed()
  p <- rasterize(p, dpi = 300)
  p
})
legend_tfs <- purrr::map(umaps_activities, function(x) {
  p <- as_ggplot(get_legend(x))
  p
})
umaps_activities <- purrr::map(umaps_activities, function(p) {
  p + theme_nothing()
})


# scATAC-seq
mat_atac <- as.matrix(df_atac[, 2:ncol(df_atac)])
rownames(mat_atac) <- df_atac$X1
colnames(mat_atac) <- c("CM", "Naive", "Non-Tfh", "Tfh")
mat_rna <- as.matrix(df_rna[, 2:ncol(df_rna)])
rownames(mat_rna) <- df_rna$X1
colnames(mat_rna) <- c("CM", "Naive", "Non-Tfh", "Tfh")
heatmaps_tfh <- purrr::map(list(mat_rna, mat_atac), function(mat) {
  ggplotify::as.ggplot(pheatmap(
    mat,
    scale = "row",
    annotation_names_row = FALSE,
    border_color = "white",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 5,
    fontsize_col = 6,
    angle_col = "45",
    legend = FALSE
  ))
})
legend_heatmap_tfh_rna <- pheatmap(mat_rna, scale = "row")
legend_heatmap_tfh_atac <- pheatmap(mat_atac, scale = "row")
# region.highlight_enhancer <- GRanges(
#   seqnames = "chr3",
#   ranges = IRanges(start = 187910000, end = 188001000)
# )
# region_bcl6_links <- "chr3-187367695-188000857"
# coverage_plot <- CoveragePlot(
#   seurat_atac,
#   region.highlight = region.highlight_enhancer,
#   group = "annotation_paper",
#   region = region_bcl6_links
# )


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
goi_th <- c("TBX21", "CXCR3", "IFNG", "RORC", "IL17F", "IL22")
umaps_th_markers <- plot_density(
  seurat_th,
  goi_th,
  reduction = "umap",
  size = 0.1
) &
  scale_color_viridis_c(option = "magma") &
  NoLegend() &
  theme(
    plot.title = element_text(size = 5, hjust = 0.5),
    axis.text = element_text(size = 0),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )


# Spatial
goi_spatial_sub <- as.character(unique(unlist(goi_spatial)))
rows_magic <- rownames(seurat_spatial[["MAGIC_Spatial"]]@data)
goi_spatial_sub <- goi_spatial_sub[goi_spatial_sub %in% rows_magic]
cor_mtrx_gg <- correlation_heatmap2( 
  se = seurat_spatial,
  genes = goi_spatial_sub,
  assay = "MAGIC_Spatial",
  slot = "data",
  text_size = 5,
  type = "upper"
) +
  scale_x_discrete(position = "top") +
  theme(
    axis.text.x = element_text(hjust = 0),
    panel.grid = element_blank()
  )
legend_cor_mtrx <- as_ggplot(get_legend(cor_mtrx_gg))
cor_mtrx_gg <- cor_mtrx_gg + NoLegend()
trajectory_gg <- SPATA2::plotTrajectory(
  object = spata_obj, 
  trajectory_name = "Tfh-migration",
  color_by = "annotation",
  pt_clrp = "npg",
  pt_size = 0.25,
  sgmt_size = 0.5,
  pt_alpha = 0.5, # reduce alpha to highlight the trajectory's course
  display_image = FALSE
) +
  SPATA2::legendTop() +
  ggplot2::scale_y_reverse() +
  theme(legend.title = element_blank())
legend_trajectory <- as_ggplot(get_legend(trajectory_gg))
trajectory_gg <- trajectory_gg + theme(legend.position = "none")
gene_vec <- as.character(unique(unlist(goi_spatial)))
gene_vec <- gene_vec[gene_vec %in% rownames(seurat_spatial)]
heat_spatial_gg <- SPATA2::plotTrajectoryHeatmap(
  object = spata_obj,
  trajectory_name = "Tfh-migration",
  variables = goi_spatial,
  arrange_rows = "maxima",
  colors = hm_colors,
  show_rownames = TRUE,
  split_columns = TRUE, 
  smooth_span = 0.5,
  border_color = NA,
  legend = FALSE,
  fontsize_row = 5
) %>%
  ggplotify::as.ggplot()
legend_heat_spatial <- as_ggplot(get_legend(heat_spatial_gg))
heat_spatial_gg <- heat_spatial_gg + theme(legend.position = "none")


# TCR
seurat_cite@reductions$wnn.umap@cell.embeddings[, 1] <- seurat_cite$CD4_UMAP1
seurat_cite@reductions$wnn.umap@cell.embeddings[, 2] <- seurat_cite$CD4_UMAP2
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
    CD4_annotation_level_5,
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


# Arrange figure
umaps <- umap_rna / (umaps_activities[[1]] | umaps_activities[[2]]) / empty_plot +
  plot_layout(heights = c(0.65, 0.35, 0.025))
top <- umaps | heatmap_rna | dot_plot_cite
top  <- top +
  plot_layout(widths = c(0.35, 0.35, 0.3))
heatmaps_tfh_arr <- heatmaps_tfh[[1]] | heatmaps_tfh[[2]]
umaps_th_arr <- umap_th / umaps_th_markers
# mid <- empty_plot
bottom1 <- cor_mtrx_gg +
  inset_element(trajectory_gg, right = 1.5, left = 0, top = 0.5, bottom = 0)
bottom <- (bottom1 | heat_spatial_gg | empty_plot) +
  plot_layout(heights = c(0.5, 0.2, 0.4))
fig <- top / bottom
tcr <- umap_tcr / barplot_tcr_gg


# Save
ggsave(
  filename = path_to_save,
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 20, 
  units = "cm"
)
ggsave(
  filename = "~/Desktop/PhD/fig2/heatmaps_tfh_fig2.pdf",
  plot = heatmaps_tfh_arr,
  device = cairo_pdf,
  width = 7, 
  height = 10,
  units = "cm"
)
ggsave(
  filename = "~/Desktop/PhD/fig2/umaps_th_arr.pdf",
  plot = umaps_th_arr,
  device = cairo_pdf,
  width = 5, 
  height = 9,
  units = "cm"
)
ggsave(
  filename = "~/Desktop/PhD/fig2/tcr_arr.pdf",
  plot = tcr,
  device = cairo_pdf,
  width = 5, 
  height = 9,
  units = "cm"
)
legends_list <- list(legend_umap, legend_heatmap, legend_dot_plot,
                     legend_tfs[[1]], legend_tfs[[2]])

names(legends_list) <- c("umap_all", "heatmap_rna", "dotplot_cite", "BCL6",
                         "PRDM1")
walk2(legends_list, names(legends_list), function(p, x) {
  ggsave(
    filename = glue::glue("~/Desktop/PhD/fig2/legends/{x}_legend_fig2_top_rows.pdf"),
    plot = p,
    device = cairo_pdf,
    width = 21, 
    height = 7, 
    units = "cm"
  )
})
