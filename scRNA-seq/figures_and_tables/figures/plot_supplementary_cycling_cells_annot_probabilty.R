# This script plots the supplementary figure containing:
# 1. UMAP (all cells) with S.Score
# 2. UMAP (all cells) with G2M Score
# 3. UMAP (all cells) with annotation probability (CITE-seq)
# 4. UMAP (all cells) with annotation probability (scATAC-seq)


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrastr)
library(here)
library(pals)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Define parameters
path_to_save <- here("results/paper/figures/supplementary_figure_cycling_cells_annotation_prob.png")
path_to_save2 <- here("results/paper/figures/supplementary_figure_umap_annotation_1.png")


# Read data
rna_multiome_cite <- read_delim(
  path_to_qc_metrics_rna_multi_cite,
  delim = ";",
  col_names = TRUE
)
atac_multiome <- read_delim(
  path_to_qc_metrics_atac_multi,
  delim = ";",
  col_names = TRUE
)
# df <- read_delim(path_to_csv, col_names = TRUE, delim = ";")


# UMAP annotation level 1
umap_annot_df <- rna_multiome_cite %>%
  filter(assay != "CITE-seq")
levels_annot_1 <- c("NBC_MBC", "GCBC", "PC", "CD4_T", "Cytotoxic",
                    "myeloid", "FDC", "epithelial", "PDC")
umap_annot_df$annotation_level_1[umap_annot_df$annotation_level_1 %in% c("preBC", "preTC")] <- "PDC"
umap_annot_df$annotation_level_1 <- factor(
  umap_annot_df$annotation_level_1,
  levels_annot_1
)
cols_annot_1 <- polychrome(length(levels_annot_1))
names(cols_annot_1) <- NULL
umap_annot_1 <- umap_annot_df %>%
  ggplot(aes(UMAP_1_level_1, UMAP_2_level_1, color = annotation_level_1)) +
  geom_point(alpha = 0.75, shape = ".") +
  scale_color_manual(values = cols_annot_1) +
  labs(title = "Annotation level 1") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 7, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) +
    guides(color = guide_legend(override.aes = list(size = 2, shape = 16))) +
    coord_fixed()
# umap_annot_1 <- rasterize(umap_annot_1, dpi = 300)


# UMAP cycling cells
(umap_s_score <-  rna_multiome_cite %>%
  ggplot(aes(UMAP_1_level_1, UMAP_2_level_1, color = S.Score)) +
    geom_point(alpha = 0.75, shape = ".") +
    scale_color_viridis_c(option = "magma") +
    labs(title = "S Score") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      legend.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank()
    ))
# umap_s_score <- rasterize(umap_s_score, dpi = 300)
(umap_g2m_score <-  rna_multiome_cite %>%
    ggplot(aes(UMAP_1_level_1, UMAP_2_level_1, color = G2M.Score)) +
    geom_point(alpha = 0.75, shape = ".") +
    scale_color_viridis_c(option = "magma") +
    labs(x = "UMAP1", y = "UMAP2", title = "G2M Score") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      legend.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank()
    ))
# umap_g2m_score <- rasterize(umap_g2m_score, dpi = 300)


# UMAPs annotation probability
umap_cite <- rna_multiome_cite %>%
  dplyr::filter(assay == "CITE-seq") %>%
  ggplot(aes(UMAP_1_level_1, UMAP_2_level_1, color = annotation_prob)) +
    geom_point(alpha = 0.9, shape = ".") +
    scale_color_distiller(palette = "Blues", direction = 1, limits = c(0, 1)) +
    labs(title = "CITE-seq", color = "P(annotation)") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
# umap_cite <- rasterize(umap_cite, dpi = 300)
umap_atac <- atac_multiome %>%
  dplyr::filter(assay == "scATAC-seq") %>%
  ggplot(aes(UMAP_1_level_1, UMAP_2_level_1, color = annotation_prob)) +
    geom_point(alpha = 0.9, shape = ".") +
    scale_color_distiller(palette = "Blues", direction = 1, limits = c(0, 1)) +
    labs(title = "scATAC-seq", color = "P(annotation)") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank()
    )
# umap_atac <- rasterize(umap_atac, dpi = 300)


# Arrange
fig <- (umap_s_score | umap_g2m_score) / (umap_cite | umap_atac)
fig <- fig &
  theme(legend.position = "none")

# Save
ggsave(
  filename = path_to_save,
  plot = fig,
  width = 20,
  height = 18,
  units = "cm",
  dpi = 300
)
ggsave(
  filename = path_to_save2,
  plot = umap_annot_1,
  width = 20,
  height = 9,
  units = "cm",
  dpi = 300
)
