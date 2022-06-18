# This script plots the relevant panels for the epithelial cells


# Load packages
library(tidyverse)
library(patchwork)
library(ggpubr)
library(Seurat)
library(here)


# Source utilities
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
epithelial <- readRDS(path_to_save_epithelial)
fdc <- readRDS(path_to_save_fdc)
spatial <- readRDS(path_to_spatial)
magic_df <- readRDS(here("data/raw_data_figures/MAGIC-mtrx-FDC.rds"))


# UMAP epithelial
epithelial_levels <- c("Basal cells", "VEGFA+", "Surface epithelium",
                       "Outer surface", "Crypt", "FDCSP epithelium")
color_palette_epithelial <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                              "#66a61e", "#e6ab02")
color_df <- data.frame(
  annotation_epithelial = epithelial_levels,
  color = color_palette_epithelial
)
umap_df_epi <- Embeddings(epithelial, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = factor(
    epithelial$annotation_20220215,
    levels = epithelial_levels
  ))
umap_epi <- ggplot(umap_df_epi, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point(size = 1, alpha = 0.85) +
  scale_color_manual(values = color_palette_epithelial) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.85, 0.75),
    legend.text = element_text(size = 6)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))


# Dot Plot
goi_epithelial <- list(
  Basal_cells = c("KRT5", "KRT14", "S100A2"),
  VEGFA = c("VEGFA", "MIR205HG", "PAX1"),
  Surface_epithelium = c("KRT4", "KRT13", "KRT78", "KRT80", "MAL", "SPRR3", "TMPRSS11B", "TMPRSS2", "ACE2"),
  Outer_surface = c("LCE3A", "LCE3D", "LCE3E", "SPRR2D", "SPRR2E", "CNFN"),
  Crypt = c("KRT8", "CD63", "IL1B", "IFI27", "S100A6", "SPIB", "MARCKSL1"),
  FDCSP_epithelium = c("FDCSP", "KRTDAP", "CALML5")
)
goi_epithelial <- unlist(goi_epithelial)
names(goi_epithelial) <- NULL
epithelial$annotation_figure <- factor(
  epithelial$annotation_20220215,
  levels = epithelial_levels
)
dot_plot_epi <- DotPlot(
  epithelial,
  features = rev(goi_epithelial),
  group.by = "annotation_figure"
)
dot_plot_epi <- dot_plot_epi +
  labs(x = "", y = "") +
  scale_size_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    range = c(0.1, 4.5)
  ) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
dot_plot_epi$guides$colour$title <- "Average\nExpression"
dot_plot_epi$guides$size$title <- "Percent\nExpressed"


# UMAP and dotplot FDC
old_levels_fdc <- c("LZ FDC", "CD14+CD55+ FDC", "DZ FDC", "FRC", "MRC",
                    "unknown")
fdc$annotation <- factor(fdc$annotation_20220215, levels = old_levels_fdc)
new_levels_fdc <- c("FDC", "CD14+CD55+ FDC", "COL27A1+ FDC", "FRC", "MRC",
                    "unknown")
levels(fdc$annotation) <- new_levels_fdc
color_palette_fdc <-  c("#632c63", "#e4624e", "#92e8df", "yellow3", "limegreen",
                        "#999999")
color_df2 <- data.frame(
  annotation_fdc = new_levels_fdc,
  color = color_palette_fdc
)
umap_fdc_df <- Embeddings(fdc, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = fdc$annotation)
umap_fdc <- ggplot(umap_fdc_df, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point(size = 1, alpha = 0.85) +
  scale_color_manual(values = color_palette_fdc) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.85, 0.75),
    legend.text = element_text(size = 6)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))
goi_fdc <- list(
  "FDC" = c("FDCSP", "CLU", "CR2", "CXCL13", "VIM"),
  "CD14+CD55+ FDC" = c("CD14", "CD55"),
  "COL27A1+ FDC" = c("COL27A1"),
  "FRC" = c("PERP", "CCL20", "S100A2", "S100A9"),
  "MRC" = c("DCN", "DIO2", "LUM", "BGN", "PDGFRB", "COL1A1", "COL1A2", "COL3A1",
            "COL5A2", "COL6A3", "COL12A1", "COL14A1", "COL18A1")
)
goi_fdc <- unlist(goi_fdc)
names(goi_fdc) <- NULL
dot_plot_fdc <- DotPlot(
  fdc,
  features = rev(goi_fdc),
  group.by = "annotation"
)
dot_plot_fdc <- dot_plot_fdc +
  labs(x = "", y = "") +
  scale_size_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    range = c(0.1, 4.5)
  ) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
dot_plot_fdc$guides$colour$title <- "Average\nExpression"
dot_plot_fdc$guides$size$title <- "Percent\nExpressed"


# Arrange 
fig_epithelial <- umap_epi | dot_plot_epi
fig_epithelial  <- fig_epithelial +
  plot_layout(widths = c(0.67, 0.33))
fig_fdc <- umap_fdc | dot_plot_fdc
fig_fdc  <- fig_fdc +
  plot_layout(widths = c(0.67, 0.33))


# Plot spatial
magic_assay <- CreateAssayObject(counts = as.matrix(magic_df))
spatial <- spatial[, colnames(magic_df)]
spatial[["MAGIC_Spatial"]] <- magic_assay
DefaultAssay(spatial) <- "MAGIC_Spatial"
goi_spatial <- c("COL1A1", "COL1A2", "COL3A1", "DCN")
image_id <- "esvq52_nluss5"
spatial_gg <- SpatialFeaturePlot(
  object = spatial,
  features = goi_spatial,
  alpha = c(0, 1),
  images = image_id,
  ncol = 2
) &
  scale_fill_gradientn(colors = heat.colors(10, rev = TRUE)) &
  scale_alpha(range = c(0, 1)) &
  theme(plot.title = element_blank(), legend.position = "none")


# Save
ggsave(
  filename =  here("results/paper/figures/supplementary_figure_epithelial_umap_dotplot.pdf"),
  plot = fig_epithelial,
  device = cairo_pdf,
  width = 16, 
  height = 12, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/supplementary_figure_FDC_umap_dotplot.pdf"),
  plot = fig_fdc,
  device = cairo_pdf,
  width = 16, 
  height = 12, 
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_FDC_spatial.png"),
  plot = spatial_gg,
  dpi = 600,
  width = 8,
  height = 8,
  units = "cm"
)
