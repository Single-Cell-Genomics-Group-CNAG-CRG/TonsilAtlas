# This script plots the relevant panels for the myeloid  cells


# Load packages
library(tidyverse)
library(patchwork)
library(ggpubr)
library(Seurat)


# Read data
seurat <- readRDS("~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/myeloid/myeloid_annotated_level_5.rds")
seurat$annotation_level_4[seurat$annotation_level_4 == "Neutrophil Granulocytes"] <- "Granulocytes" 


# UMAP
myeloid_levels <- c("DC1 precursor", "DC1 mature", "DC2", "DC3", "DC4", "DC5",
                    "IL7R DC", "aDC1", "aDC2", "aDC3", "IL7R MMP12 macrophages",
                    "ITGAX ZEB2 macrophages", "C1Q HLA macrophages",
                    "SELENOP FUCA1 PTGDS macrophages", "M1 Macrophages",
                    "Monocytes", "Mast cells", "Granulocytes",
                    "Cycling")
myeloid_colors <- c("#ff8fd2", "#fc7da7", "#f6747f", "#f35449", "#f25e31",
                    "#f08819", "#e1b40e", "#e2f155", "#b4d84b", "#6ead34",
                    "#91f8ca", "#7ef1eb", "#68cde8", "#519bdb", "#3f56ca",
                    "#cd87de", "#c035a7", "black", "gray")
seurat$annotation_figure_3 <- factor(
  seurat$annotation_level_4,
  levels = myeloid_levels
)
color_df <- data.frame(
  annotation_myeloid = myeloid_levels,
  color = myeloid_colors
)
umap_df <- Embeddings(seurat, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = factor(
    seurat$annotation_figure_3,
    levels = myeloid_levels
  ))
umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point(size = 1, alpha = 0.85) +
  scale_color_manual(values = myeloid_colors) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))


# Dot Plot
goi_epithelial <- unlist(goi_epithelial)
names(goi_epithelial) <- NULL
seurat$annotation_figure <- factor(
  seurat$annotation_level_3,
  levels = epithelial_levels
)
dot_plot <- DotPlot(
  seurat,
  features = rev(goi_epithelial),
  group.by = "annotation_figure"
)
dot_plot <- dot_plot +
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
dot_plot$guides$colour$title <- "Average\nExpression"
dot_plot$guides$size$title <- "Percent\nExpressed"


# Plot signatures spatial transcriptomics: WRONG!!! Use spatial deconvolution instead
# spatial <- subset(spatial, gem_id == "tarwe1_xott6q")
# spatial@images <- spatial@images[Seurat::Images(spatial) == "tarwe1_xott6q"]
# DefaultAssay(spatial) <- "Spatial"
# spatial <- AddModuleScore(
#   spatial,
#   features = goi_epithelial[c("Basal_cells", "Surface_epithelium", "Crypt")],
#   name = c("Basal_cells", "Surface_epithelium", "Crypt")
# )
# SpatialFeaturePlot(spatial, features = c("Basal_cells1", "Surface_epithelium2", "Crypt3"))



# Arrange
fig <- umap | dot_plot
fig  <- fig +
  plot_layout(widths = c(0.67, 0.33))


# Save
ggsave(
  filename = "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/plots/paper/non_lymphoid/epithelial_umap_and_dot_plot.pdf",
  plot = fig,
  device = cairo_pdf,
  width = 16, 
  height = 12, 
  units = "cm"
)

