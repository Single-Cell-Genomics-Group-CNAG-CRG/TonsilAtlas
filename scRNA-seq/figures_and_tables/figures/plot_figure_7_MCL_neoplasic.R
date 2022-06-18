# This script plots figure related with the neoplasic cells
# of mantle cell lymphoma (MCL)


# Load packages
library(Seurat)
library(Signac)
library(ggrastr)
library(ggpubr)
library(tidyverse)
library(here)


# Read data
seurat_102 <- readRDS(here("MCL/results/R_objects/7.seurat_tumoral_102_clustered.rds"))


# UMAP annotation
colors_mcl <- c(
  "chrY+ metallothionein" = "#f56666",
  "chrY+ CD69+AP1+" = "#6df874",
  "chrY+ MIR155HG+NFKB1+MYC+" = "#87d0f7",
  "chrY+ RGS7+BANK1+" = "#c737d7",
  "chrY- metallothionein" = "#ab2121",
  "chrY- CD69+AP1+" = "#36a13b",
  "chrY- MIR155HG+NFKB1+MYC+" = "#3f82a6",
  "Cycling tumoral" = "gray",
  "non-tumoral B-cells" = "#d78437"
)
seurat_102$annotation <- factor(
  seurat_102$annotation_20220518,
  levels = names(colors_mcl)
)
umap_title <- format(ncol(seurat_102), big.mark = ",", scientific = FALSE)
umap_annotation <- Embeddings(seurat_102, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = seurat_102$annotation) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point(shape = ".", alpha = 0.9) +
  scale_color_manual(values = colors_mcl, breaks = names(colors_mcl)) +
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
    legend.text = element_text(size = 5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))
legend_umap_annotation <- as_ggplot(get_legend(umap_annotation))
umap_annotation <- umap_annotation + theme(legend.position = "none")
umap_annotation <- ggrastr::rasterize(umap_annotation, dpi = 600)


# UMAP chrY
cols_chrY <- c(
  "chrY-" = "#d03132",
  "chrY+" = "azure4"
)
umap_chrY <- Embeddings(seurat_102, "umap") %>%
  as.data.frame() %>%
  mutate(has_loss_chrY = seurat_102$has_loss_chrY) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = has_loss_chrY)) +
  geom_point(shape = ".") +
  scale_color_manual(values = cols_chrY, breaks = names(cols_chrY)) +
  labs(x = "UMAP1", y = "UMAP2") +
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
    legend.text = element_text(size = 6),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))
legend_umap_chrY <- as_ggplot(get_legend(umap_chrY))
umap_chrY <- umap_chrY + theme(legend.position = "none")
umap_chrY <- ggrastr::rasterize(umap_chrY, dpi = 600)


# Plot dot plot
goi_102 <- c("SOX11", "CCND1",
             "MT2A", "MT1X", "MT1G",
             "JUN", "FOS", "CD69", "RGS7", "BANK1",
             "MIR155HG", "TRAF1", "NFKB1", "NFAT5", "MYC", "IL2RA",
             "TOP2A")
seurat_102$annotation2 <- factor(
  seurat_102$annotation_20220518,
  rev(names(colors_mcl))
)
Idents(seurat_102) <- "annotation2"
dot_plot <- DotPlot(seurat_102, features = goi_102) +
  # coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  scale_size_continuous(range = c(0.1, 3.5)) +
  theme(
    legend.position = "right",
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    axis.text.x =  element_text(size = 6, angle = 90, vjust = 1, hjust = 1),
    # axis.text = element_blank()
  )
# legend_dot_plot <- as_ggplot(get_legend(dot_plot))
# dot_plot <- dot_plot + theme(legend.position = "none")

# Plot violin plot cycling cells
s_score_gg <- seurat_102@meta.data %>%
  ggplot(aes(annotation, S.Score, color = has_loss_chrY)) +
   geom_jitter(size = 0.01) +
   scale_color_manual(values = cols_chrY, breaks = names(cols_chrY)) +
   theme_classic() +
   theme(
     axis.text.x = element_blank(),
     axis.text.y = element_text(size = 6),
     axis.title.x = element_blank(),
     axis.title.y = element_text(size = 6),
     legend.position = "none")
s_score_gg <- rasterize(s_score_gg, dpi = 600)
g2m_score_gg <- seurat_102@meta.data %>%
  ggplot(aes(annotation, G2M.Score, color = has_loss_chrY)) +
  geom_jitter(size = 0.01) +
  scale_color_manual(values = cols_chrY, breaks = names(cols_chrY)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    legend.position = "none")
g2m_score_gg <- rasterize(g2m_score_gg, dpi = 600)
cc_gg <- s_score_gg / g2m_score_gg

# Plot heatmap de novo expression in MCL



# Save
ggsave(
  filename = here("results/paper/figures/figure_7_MCL_umap_neoplasic_annotation.pdf"),
  plot = umap_annotation,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_7_MCL_umap_neoplasic_chrY.pdf"),
  plot = umap_chrY,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_7_MCL_umap_neoplasic_annotation_legend.pdf"),
  plot = legend_umap_annotation,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_7_MCL_umap_neoplasic_chrY_legend.pdf"),
  plot = legend_umap_chrY,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_7_MCL_dotplot_neoplasic.pdf"),
  plot = dot_plot,
  width = 12,
  height = 6,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/figure_7_MCL_jitter_cell_cycle.pdf"),
  plot = cc_gg,
  width = 5,
  height = 5.5,
  units = "cm"
)
