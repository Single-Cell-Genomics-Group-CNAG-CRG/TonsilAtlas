# This script plots supplementary figure related with the neoplasic cells
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
seurat_413 <- readRDS(here("MCL/results/R_objects/7.seurat_tumoral_413_clustered.rds"))


# UMAPs genes expressed in chrY+
goi_chrY <- c("UTY", "KDM5D", "DDX3Y", "USP9Y", "ZFY", "EIF1AY")
umaps_chrY <- purrr::map(list(seurat_102, seurat_413), function(seurat_obj) {
  plots <- purrr::map(goi_chrY, function(x) {
    p <- Embeddings(seurat_obj, "umap") %>%
      as.data.frame() %>%
      mutate(expression = seurat_obj[["RNA"]]@data[x, ]) %>% 
      ggplot(aes(UMAP_1, UMAP_2, color = expression)) +
      geom_point(shape = ".") +
      scale_color_distiller(palette = "Blues", direction = 1) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 6, face = "plain"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) +
      coord_fixed()
    p
  })
  plots
})
legend_umaps_chrY <- as_ggplot(get_legend(umaps_chrY[[1]][[1]]))
rm(seurat_102)
gc()
umaps_chrY <- c(umaps_chrY[[1]], umaps_chrY[[2]])
umaps_chrY <- purrr::map(umaps_chrY, function(p) {
  p + theme(legend.position = "none")
})
umaps_chrY_arranged <- ggarrange(plotlist = umaps_chrY, nrow = 2, ncol = 6)


# Density plots + chry Score
density_102 <- seurat_102@meta.data %>%
  ggplot(aes(chrY_score1)) +
  geom_density() +
  xlab("chrY expression Score") +
  geom_vline(xintercept = -0.32, linetype = "dashed", color = "darkblue") +
  theme_classic() +
  theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6))
density_413 <- seurat_413@meta.data %>%
   ggplot(aes(chrY_score1)) +
     geom_density() +
     xlab("chrY expression Score") +
     geom_vline(xintercept = 0, linetype = "dashed", color = "darkblue") +
     theme_classic() +
  theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6))
umap_chrYscores <- purrr::map(list(seurat_102, seurat_413), function(seurat_obj) {
  p <- Embeddings(seurat_obj, "umap") %>%
    as.data.frame() %>%
    mutate(chrY_score = seurat_obj$chrY_score1) %>% 
    ggplot(aes(UMAP_1, UMAP_2, color = chrY_score)) +
    geom_point(shape = ".") +
    scale_color_distiller(palette = "Blues", direction = 1) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 6, face = "plain"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none") +
    coord_fixed()
  rasterize(p, dpi = 300)
})
chrY_score_fig <-
  density_102 | umap_chrYscores[[1]] | density_413 | umap_chrYscores[[2]]



# UMAP annotation
colors_mcl <- c(
  "chrY+ TSHZ2+MARCH1+" = "#f56666",
  "chrY+ CD69+JUNB+CXCR4+" = "#6df874",
  "chrY+ MIR155HG+NFKB1+MYC+" = "#87d0f7",
  "chrY- PCDH9+MARCH1+" = "#ab2121",
  "chrY- CD69+JUNB+CXCR4+" = "#36a13b",
  "chrY- MIR155HG+NFKB1+MYC+" = "#3f82a6",
  "Cycling tumoral" = "gray",
  "non-tumoral B-cells" = "#d78437",
  "undefined" = "beige"
)
seurat_413$annotation <- factor(
  seurat_413$annotation_20220523,
  levels = names(colors_mcl)
)
umap_title <- format(ncol(seurat_413), big.mark = ",", scientific = FALSE)
umap_annotation <- Embeddings(seurat_413, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = seurat_413$annotation) %>% 
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
umap_chrY <- Embeddings(seurat_413, "umap") %>%
  as.data.frame() %>%
  mutate(has_loss_chrY = seurat_413$has_loss_chrY) %>% 
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
goi_413 <- c("SOX11", "CCND1", "TSHZ2", "CCSER1", "MARCH1", "PRDM2", "CXCR4", "CD83",
             "JUNB", "CD69", "MIR155HG", "TRAF1", "NFKB1", "NFAT5",
             "MYC", "IL2RA", "PCDH9", "CCND3", "CD38", "TOP2A")
seurat_413$annotation2 <- factor(
  seurat_413$annotation,
  rev(names(colors_mcl))
)
Idents(seurat_413) <- "annotation2"
dot_plot <- DotPlot(seurat_413, features = goi_413) +
  # coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  scale_size_continuous(
    limits = c(0, 100),
    range = c(0.01, 4),
    breaks = c(0, 25, 50, 75, 100)
  ) +
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
s_score_gg <- seurat_413@meta.data %>%
  ggplot(aes(annotation, S.Score, color = has_loss_chrY)) +
  geom_jitter(size = 0.25) +
  scale_color_manual(values = cols_chrY, breaks = names(cols_chrY)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    legend.position = "none")
s_score_gg <- rasterize(s_score_gg, dpi = 300)
g2m_score_gg <- seurat_413@meta.data %>%
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
g2m_score_gg <- rasterize(g2m_score_gg, dpi = 300)
cc_gg <- s_score_gg / g2m_score_gg

# Plot heatmap de novo expression in MCL



# Save
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_umap_chrY.png"),
  plot = umaps_chrY_arranged,
  dpi = 300,
  # device = "cairo",
  width = 20,
  height = 7,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_umap_chrY_legend.pdf"),
  plot = legend_umaps_chrY,
  width = 20,
  height = 7,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_chrY_scores.pdf"),
  plot = chrY_score_fig,
  width = 20,
  height = 7,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_umap_neoplasic_annotation_413.pdf"),
  plot = umap_annotation,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_umap_neoplasic_chrY.pdf"),
  plot = umap_chrY,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_umap_neoplasic_413_annotation_legend.pdf"),
  plot = legend_umap_annotation,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_umap_neoplasic_chrY_legend.pdf"),
  plot = legend_umap_chrY,
  width = 5,
  height = 5,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_dotplot_neoplasic.pdf"),
  plot = dot_plot,
  width = 13.5,
  height = 7,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_MCL_jitter_cell_cycle.pdf"),
  plot = cc_gg,
  width = 5,
  height = 5.5,
  units = "cm"
)
