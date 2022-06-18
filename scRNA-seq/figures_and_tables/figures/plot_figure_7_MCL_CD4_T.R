# This script plots figure related with the TME of mantle cell lymphoma (MCL, CD4 T)


# Load packages
library(Seurat)
library(ggrastr)
library(tidyverse)
library(here)


# Source utilities
source(here("scRNA-seq/bin/utils_figure2.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
path_to_cd4_t_102 <- here("MCL/results/R_objects/7.seurat_CD4_T_102_annotated.rds")
path_to_cd4_t_413 <- here("MCL/results/R_objects/7.seurat_CD4_T_413_annotated.rds")
cd4_tonsil <- readRDS(path_to_save_cd4)
cd4_102 <- readRDS(path_to_cd4_t_102)
cd4_413 <- readRDS(path_to_cd4_t_413)
dea_tregs <- readRDS("data/raw_data_figures/dea_tregs_list.rds")


# Plot UMAPs age groups
cd4_tonsil$annotation_paper <- factor(
  cd4_tonsil$annotation_20220215,
  levels = names(colors_rna)
)
age_groups2 <- list("kid", c("young_adult", "old_adult"))
n_cells <- table(cd4_tonsil$age_group)
n_cells2 <- c(
  "kid" = n_cells["kid"],
  "adult" = (n_cells["young_adult"] + n_cells["old_adult"])
)
umaps_age_groups_cd4 <- purrr::map2(age_groups2, n_cells2, function(x, y) {
  print(x)
  p_title <- format(y, big.mark = ",", scientific = FALSE)
  p <- Embeddings(cd4_tonsil, "umap") %>%
    as.data.frame() %>%
    dplyr::mutate(
      annotation_paper = cd4_tonsil$annotation_paper,
      age_group = cd4_tonsil$age_group
    ) %>%
    dplyr::filter(age_group %in% x) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = annotation_paper)) +
    geom_point(shape = ".", alpha = 0.85) +
      scale_color_manual(values = colors_rna, breaks = names(colors_rna)) +
      theme_classic() +
      labs(
        title = str_c(x, collapse = " "),
        subtitle = str_c(p_title, "cells", sep = " ")
      ) +
      coord_fixed() +
      theme(
        # legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 7),
        plot.subtitle = element_text(size = 6),
        axis.title = element_blank(),
        # legend.text = element_text(size = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      ) +
      guides(
        color = guide_legend(
          override.aes = list(shape = 16, size = 2),
          nrow = 2,
          byrow = FALSE
        )
      )
    p <- rasterize(p, dpi = 300)
  p
})


# Plot UMAPs MCL
cd4_102$annotation_paper <- factor(
  cd4_102$annotation_SLOcatoR,
  levels = names(colors_rna)
)
cd4_413$annotation_paper <- factor(
  cd4_413$annotation_SLOcatoR,
  levels = names(colors_rna)
)
mcl_metadata <- rbind(cd4_102@meta.data, cd4_413@meta.data)
p_title <- format(nrow(mcl_metadata), big.mark = ",", scientific = FALSE)
umap_mcl_cd4 <- mcl_metadata %>%
  ggplot(aes(UMAP_1_SLOcatoR, UMAP_2_SLOcatoR, color = annotation_paper)) +
  geom_point(shape = ".", alpha = 1) +
  scale_color_manual(values = colors_rna, breaks = names(colors_rna)) +
  theme_classic() +
  labs(
    title = "MCL1 and MCL2",
    subtitle = str_c(p_title, "cells", sep = " ")
  ) +
  coord_fixed() +
  theme(
    # legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 7),
    plot.subtitle = element_text(size = 6),
    axis.title = element_blank(),
    # legend.text = element_text(size = 6),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 16, size = 2),
      nrow = 2,
      byrow = FALSE
    )
  )
umap_mcl_cd4 <- rasterize(umap_mcl_cd4, dpi = 300)


# Barplot
cd4_tonsil$group <- cd4_tonsil$age_group
cd4_102$group <- "MCL1"
cd4_413$group <- "MCL2"
barplot_dfs <- purrr::map(list(cd4_tonsil, cd4_102, cd4_413), function(seurat_obj) {
  df <- seurat_obj@meta.data[, c("annotation_paper", "group")]
  df
})
new_levels <- c("kid", "young_adult", "old_adult", "MCL1", "MCL2")
barplot_df <- barplot_dfs %>%
  bind_rows() %>%
  group_by(group, annotation_paper) %>%
  summarise(n_cells = n()) %>%
  mutate(
    n_cells_total = sum(n_cells),
    percentage_cells = round(n_cells / n_cells_total * 100, 3),
    group = factor(group, levels = new_levels)
  )
barplot_gg <- ggplot(barplot_df, aes(group, percentage_cells, fill = annotation_paper)) +
  geom_col() +
  scale_fill_manual(values = colors_rna, breaks = names(colors_rna)) +
  ylab("Percentage of cells (%)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 6, color = "black")
  )
barplot_gg


# Volcano plot
top_genes <- c("MT2A", "MT1X", "MT1E", "STAT4", "TNFRSF9", "CXCR4", "IKZF2",
               "TSHZ2", "LAG3", "SLA", "IL7R", "MAF", "IL1R1")
dea_sub <- dea_tregs$`Eff-Tregs`[top_genes, ]
volcano_gg <- dea_tregs$`Eff-Tregs` %>%
    ggplot(aes(avg_log2FC, -1 * log10(p_val_adj), color = is_significant)) +
    geom_point(size = 0.25) +
    geom_text_repel(data = dea_sub, aes(label = gene), color = "black",
                    size = 1.75, max.overlaps = 15) +
    scale_color_manual(
      values = c("gray", "forestgreen"),
      breaks = c("not significant", "significant"),
      labels = c("not sig.", "sig.")
    ) +
    labs(x =  "log2 (MCL / HD)", y = "-log10(adj. p-value)") +
    theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "top"
    ) +
    guides(color = guide_legend(override.aes = list(size = 2)))


# Arrange
umaps <- umaps_age_groups_cd4[[1]] | umaps_age_groups_cd4[[2]] | umap_mcl_cd4


# Save
ggsave(
  plot = umaps,
  filename = here("results/paper/figures/figure_7_MCL_umaps_CD4_T.pdf"),
  width = 13,
  height = 8,
  units = "cm"
)
ggsave(
  plot = barplot_gg,
  filename = here("results/paper/figures/figure_7_MCL_barplot_CD4_T.pdf"),
  width = 4,
  height = 5,
  units = "cm"
)
ggsave(
  plot = volcano_gg,
  filename = here("results/paper/figures/figure_7_MCL_volcano_plot_CD4_T.pdf"),
  width = 4.25,
  height = 6,
  units = "cm"
)

