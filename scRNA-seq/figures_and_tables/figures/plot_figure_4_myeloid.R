# This script plots all panels related with myeloid cells for figure 4

# Load packages
library(Seurat)
library(scCustomize)
library(pheatmap2)
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(here)
library(ggpubr)
library(ggrastr)
library(textshape)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure4.R"))


# Read data
seurat <- readRDS(path_to_save_myeloid)
spatial <- readRDS(path_to_spatial)
magic_df <- readRDS(here("data/raw_data_figures/MAGIC-mtrx.rds"))
gsea_list <- readRDS(here("scRNA-seq/results/R_objects/gsea_list_slancytes.rds"))
slan_markers <- read.delim(here("data/raw_data_figures/BA2020Faseb_SLANmarkers.txt"))


# Plot UMAP
seurat$annotation_paper <- factor(
  seurat$annotation_20220215,
  levels = names(colors_myel)
)
Idents(seurat) <- "annotation_paper"
umap_title <- format(ncol(seurat), big.mark = ",", scientific = FALSE)
umap_annotation <- seurat@meta.data %>%
  ggplot(aes(UMAP_1_20220215, UMAP_2_20220215, color = annotation_paper)) +
  geom_point(size = 0.25) +
  scale_color_manual(values = colors_myel, breaks = names(colors_myel)) +
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
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))
legend_umap <- as_ggplot(get_legend(umap_annotation))
umap_annotation <- umap_annotation + theme(legend.position = "none")
umap_annotation <- ggrastr::rasterize(umap_annotation, dpi = 600)


# Markers
dot_plot <- DotPlot(seurat, features = rev(goi_myeloid)) +
  coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    axis.text.x = element_blank())
legend_dot_plot <- as_ggplot(get_legend(dot_plot))
dot_plot <- dot_plot + theme(legend.position = "none")


# Markers slancytes
slan_levels <- c("IL7R MMP12 macrophages", "C1Q HLA macrophages",
                 "SELENOP FUCA1 PTGDS macrophages", "ITGAX ZEB2 macrophages")
seurat_slan <- subset(seurat, idents = slan_levels)
seurat_slan$annotation_paper <- factor(
  seurat_slan$annotation_paper,
  rev(slan_levels)
)
Idents(seurat_slan) <- "annotation_paper"
goi_slan <- c("SLAMF7", "FCGR2A", "TNFRSF1B", "IL17RA", "SIGLEC10", "ICAM1", "LILRB1",
              "LILRB3", "MMP9", "MMP12", "MMP14", "MMP25", "TLR1", "TLR2",
              "STAB1", "SCARB1", "SCARB2",
              "C1QA", "C1QB", "C1QC", "HLA-DRA",
              "HLA-DPB1", "HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "S100B",
              "SELENOP", "FOLR2", "FUCA1", "PTGDS", "APOC1", "APOE",
              "CCL18", "ATP5F1E", "UQCRB",
              "ITGAX", "ZEB2", "NR1H3", "IL6R", "IL4R", "IL18", "JAK3", "STAT2", "CR1",
              "NFATC2IP", "CD163L1",
              "MAF", "MAFB", "MAFG")

(dot_plot_slan <- DotPlot(seurat_slan, features = goi_slan) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  scale_size_continuous(range = c(0.1, 3.5)) +
  theme(
    axis.title = element_blank(),
    # axis.text.y = element_text(size = 6),
    axis.text.y = element_blank(),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 5, angle= 45, hjust = 1, vjust = 1)))
legend_dotplot_slan <- as_ggplot(get_legend(dot_plot_slan))
dot_plot_slan <- dot_plot_slan + theme(legend.position = "none")

# GSEA plots
gsea_list <- gsea_list[slan_levels]
sel_terms <- c(
  "oxidative phosphorylation" = "GO:0006119",
  "MHC class II protein complex assembly" = "GO:0002399",
  "extracellular matrix disassembly" = "GO:0022617",
  "complement activation, classical pathway" = "GO:0006958",
  "tumor necrosis factor superfamily cytokine production "= "GO:0071706"
)
gsea_plots <- purrr::map(names(sel_terms), function(x) {
  plots <- purrr::map(gsea_list, function(obj) {
    p <- gseaplot(obj, by = "runningScore", geneSetID = sel_terms[x])
    p +
      scale_y_continuous(limits = c(-1, 1)) +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
  })
  patchwork::wrap_plots(plots, ncol = 4)
})
gsea_arranged <- patchwork::wrap_plots(gsea_plots, nrow = 5)


# Heatmap slan markers Bianchetto-Aguilera
new_order <- c("IL7R MMP12 macrophages", "C1Q HLA macrophages",
               "SELENOP FUCA1 PTGDS macrophages", "ITGAX ZEB2 macrophages",
               "DC1 precursor", "DC1 mature", "DC2", "DC3", "DC4", "DC5", "IL7R DC",
               "aDC1", "aDC2", "aDC3", "M1 Macrophages", "Monocytes",
               "Mast cells", "Neutrophil Granulocytes", "Cycling")
colors_myel <- colors_myel[new_order]
sel_genes <- slan_markers$Gene[slan_markers$Gene %in% rownames(seurat)]
avgexpr_mat <- AverageExpression(
  seurat,
  features = sel_genes,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "annotation_paper",
  slot = "data"
)$RNA
avgexpr_mat <- avgexpr_mat[, new_order]
mycolors <- list(cell_type = colors_myel)
annotation_col <- data.frame(cell_type = names(colors_myel)) 
rownames(annotation_col) <- names(colors_myel)
zero_rows <- map_lgl(1:nrow(avgexpr_mat), \(.) diff(range(avgexpr_mat[., ])) == 0)
avgexpr_mat <- avgexpr_mat[!zero_rows, ]
input_mat <- t(apply(avgexpr_mat, 1, function(x) (x - min(x)) / diff(range(x))))
tmp <- slan_markers$Gene[slan_markers$CellType == "SLAN"]
slan_avg_mat <- input_mat[rownames(input_mat) %in% tmp, slan_levels]
slan_rows <- rownames(textshape::cluster_matrix(slan_avg_mat, dim = "row"))
dc_levels <- c("aDC1", "aDC2", "aDC3", "DC1 mature", "DC1 precursor", "DC2",
               "DC3", "DC4", "DC5", "IL7R DC")
tmp <- slan_markers$Gene[slan_markers$CellType == "DC"]
dc_avg_mat <- input_mat[rownames(input_mat) %in% tmp, dc_levels]
dc_rows <- rownames(textshape::cluster_matrix(dc_avg_mat, dim = "row"))
tmp <- slan_markers$Gene[slan_markers$CellType == "MAC"]
mac_avg_mat <- input_mat[rownames(input_mat) %in% tmp, "M1 Macrophages", drop = FALSE]
mac_rows <- names(textshape::cluster_matrix(mac_avg_mat, dim = "row"))
input_mat <- input_mat[c(slan_rows, dc_rows, mac_rows), ]
in2mm <- 25.4
path_save_heatmap_slan <- here("results/paper/figures/figure_4_heatmap_myeloid_slan.pdf")
pdf(path_save_heatmap_slan, width = (55 / in2mm), height = (45 / in2mm), paper = "special")
pheatmap2(
  input_mat,
  annotation_col = annotation_col,
  annotation_colors = mycolors,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE, 
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = c(length(slan_rows), length(slan_rows) + length(dc_rows)),
  gaps_col = c(4, 14, 15),
  fontsize_row = 6
)
dev.off()


# Plot spatial
magic_assay <- CreateAssayObject(counts = as.matrix(magic_df))
spatial <- spatial[, colnames(magic_df)]
spatial[["MAGIC_Spatial"]] <- magic_assay
DefaultAssay(spatial) <- "MAGIC_Spatial"
goi_spatial <- c("MMP12", "C1QA", "SELENOP")
image_id <- "esvq52_nluss5"
spatial_gg <- SpatialFeaturePlot(
  object = spatial,
  features = goi_spatial,
  alpha = c(0, 1),
  images = image_id
) &
  scale_fill_gradientn(colors = heat.colors(10, rev = TRUE)) &
  scale_alpha(range = c(0, 1)) &
  theme(plot.title = element_blank())

# Save
ggsave(
  filename =  here("results/paper/figures/figure_4_umap_myeloid.pdf"),
  plot = umap_annotation,
  width = 7, 
  height = 8.5, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/figure_4_umap_myeloid_legend.pdf"),
  plot = legend_umap,
  width = 12, 
  height = 8.5, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/figure_4_dotplot_markers.pdf"),
  plot = dot_plot,
  width = 12.5, 
  height = 8, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/figure_4_dotplot_markers_legend.pdf"),
  plot = legend_dot_plot,
  width = 12.5, 
  height = 8, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/figure_4_dotplot_markers_slan.pdf"),
  plot = dot_plot_slan,
  width = 18, 
  height = 5, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/figure_4_dotplot_markers_slan_legend.pdf"),
  plot = legend_dotplot_slan,
  width = 18, 
  height = 5, 
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_4_gsea.pdf"),
  plot = gsea_arranged,
  width = 16,
  height = 20,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_4_spatial.png"),
  plot = spatial_gg,
  dpi = 600,
  width = 21,
  height = 7,
  units = "cm"
)
