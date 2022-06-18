# This script plots the main UMAPS of figure 1


# Load packages
library(tidyverse)
library(patchwork)
library(here)
library(ggrastr)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure1.R"))


# Read data
df <- read_delim(path_to_df_umap_fig1, delim = ";", col_names = TRUE)
df <- df %>%
  dplyr::filter(annotation_prob > 0.6 | is.na(annotation_prob))


# UMAP annotation figure 1
umap_title <- format(nrow(df), big.mark = ",", scientific = FALSE)
umap_annotation <- df %>%
  ggplot(aes(UMAP_1_level_1, UMAP_2_level_1, color = annotation_figure_1)) +
    geom_point(shape = ".", alpha = 0.9) +
    scale_color_manual(values = cols_fig1, breaks = names(cols_fig1)) +
    labs(title = str_c(umap_title, "cells", sep = " "), x = "UMAP1", y = "UMAP2") +
    theme_classic() +
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
    guides(color = guide_legend(override.aes = list(shape = 16, size = 2)))
umap_annotation <- rasterize(umap_annotation, dpi = 300)


# UMAP split by assay
n_cells <- purrr::map_dbl(cols_assays_df$assay, function(x) {
  out <- sum(df$assay == x)
  out
})
names(n_cells) <- cols_assays_df$assay
umaps_list <- purrr::map(cols_assays_df$assay, function(x) {
  df2 <- df[df$assay == x, ]
  n <- format(n_cells[x], big.mark = ",", scientific = FALSE)
  p <- df2 %>%
    ggplot(aes(UMAP_1_level_1, UMAP_2_level_1)) +
    geom_point(shape = ".", color = cols_assays_df[cols_assays_df$assay == x, "color"]) +
    labs(title = x, subtitle = str_c(n, "cells", sep = " ")) +
    theme_classic() +
    labs(color = "") +
    theme_nothing2() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 7.5),
      plot.subtitle = element_text(size = 6.5),
      legend.position = "none"
    )
    coord_fixed(ratio = 1)
  rasterize(p, dpi = 300)
})
umaps_assay_arr <- (umaps_list[[1]] | umaps_list[[2]]) /
                   (umaps_list[[3]] | umaps_list[[4]])


# Arrange
fig <- umap_annotation | umaps_assay_arr


# Save
ggsave(
  plot = fig,
  filename = path_to_save_umaps_fig1,
  width = 19.5,
  height = 12,
  units = "cm"
)
