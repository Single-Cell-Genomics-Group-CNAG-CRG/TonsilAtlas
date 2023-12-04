# This script plots the results of the cell cycle correction

# Load packages
library(tidyverse)
library(ggpubr)
library(here)
library(ggrastr)


# Read, plot, save
umap_annot_original <- readRDS(here("scRNA-seq/6-revision/tmp/umap_annot_original.rds"))
umap_cc_original <- readRDS(here("scRNA-seq/6-revision/tmp/umap_cc_original.rds"))
umap_annot_original <- umap_annot_original +
  theme(
    legend.position = "top",
    legend.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(0, "cm"),
    legend.key.width = unit(0, "cm")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 2),
      nrow = 4,
      byrow = FALSE
    )
  )
legend_annotation <- as_ggplot(get_legend(umap_annot_original))
umap_annot_original <- umap_annot_original +
  theme(legend.position = "none")
umap_cc_original <- umap_cc_original +
  theme(legend.position = "top") +
  guides(
    color = guide_legend(
      override.aes = list(size = 2),
      nrow = 1,
      byrow = FALSE
    )
  )
legend_cc <- as_ggplot(get_legend(umap_cc_original))
umap_cc_original <- umap_cc_original +
  theme(legend.position = "none")
umaps_original <- umap_annot_original / umap_cc_original
ggsave(
  plot = umaps_original,
  filename = here("results/paper/figures/revision_cell_cycle_correction_original.pdf"),
  width = 7,
  height = 14,
  units = "cm"
)
ggsave(
  plot = legend_annotation,
  filename = here("results/paper/figures/revision_cell_cycle_correction_legend_annotation.pdf"),
  width = 20,
  height = 7,
  units = "cm"
)
ggsave(
  plot = legend_cc,
  filename = here("results/paper/figures/revision_cell_cycle_correction_legend_cc.pdf"),
  width = 20,
  height = 7,
  units = "cm"
)
rm(umaps_original, umap_annot_original, umap_cc_original)
gc()

umap_annot_no_cc_genes <- readRDS(here("scRNA-seq/6-revision/tmp/umap_annot_no_cc_genes.rds"))
umap_cc_no_cc_genes <- readRDS(here("scRNA-seq/6-revision/tmp/umap_cc_no_cc_genes.rds"))
umap_annot_no_cc_genes <- umap_annot_no_cc_genes +
  theme(legend.position = "none")
umap_cc_no_cc_genes <- umap_cc_no_cc_genes +
  theme(legend.position = "none")
umap_annot_no_cc_genes <- rasterise(umap_annot_no_cc_genes, dpi = 300)
umap_cc_no_cc_genes <- rasterise(umap_cc_no_cc_genes, dpi = 300)
umaps_no_cc_genes <- umap_annot_no_cc_genes / umap_cc_no_cc_genes
ggsave(
  plot = umaps_no_cc_genes,
  filename = here("results/paper/figures/revision_cell_cycle_correction_no_cc_genes.pdf"),
  width = 7,
  height = 14,
  units = "cm"
)
rm(umap_annot_no_cc_genes, umap_cc_no_cc_genes, umaps_no_cc_genes)
gc()

umap_annot_regressed <- readRDS(here("scRNA-seq/6-revision/tmp/umap_annot_regressed.rds"))
umap_annot_regressed <- umap_annot_regressed +
  theme(legend.position = "none")
umap_annot_regressed <- rasterise(umap_annot_regressed, dpi = 300)
umap_cc_regressed <- readRDS(here("scRNA-seq/6-revision/tmp/umap_cc_regressed.rds"))
umap_cc_regressed <- umap_cc_regressed +
  theme(legend.position = "none")
umap_cc_regressed <- rasterize(umap_cc_regressed, dpi = 300)
umaps_regressed <- umap_annot_regressed / umap_cc_regressed
ggsave(
  plot = umaps_regressed,
  filename = here("results/paper/figures/revision_cell_cycle_correction_regressed.pdf"),
  width = 7,
  height = 14,
  units = "cm"
)
rm(umap_annot_regressed, umap_cc_regressed, umaps_regressed)
gc()

lisi_original <- readRDS(here("scRNA-seq/6-revision/tmp/lisi_original.rds"))
lisi_no_cc_genes <- readRDS(here("scRNA-seq/6-revision/tmp/lisi_no_cc_genes.rds"))
lisi_regressed <- readRDS(here("scRNA-seq/6-revision/tmp/lisi_regressed.rds"))
lisi_df <- bind_rows(list(lisi_original, lisi_no_cc_genes, lisi_regressed))
lisi_df$group <- factor(
  lisi_df$group,
  levels = c("No CC correction", "Removed CC genes", "CC regressed out")
)
lisi_annot_gg <- lisi_df  %>%
  ggplot(aes(group, annotation_20230124)) +
  geom_boxplot(outlier.size = 0.1) +
  ylab("LISI (Annotation)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 6))
lisi_annot_gg <- rasterise(lisi_annot_gg, dpi = 300)
lisi_cc_phase_gg <- lisi_df  %>%
  ggplot(aes(group, Phase)) +
    geom_boxplot(outlier.size = 0.1) +
    ylab("LISI (CC phase)") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 6))
lisi_cc_phase_gg <- rasterize(lisi_cc_phase_gg, dpi = 300)
lisi_gg <- lisi_annot_gg / lisi_cc_phase_gg
ggsave(
  plot = lisi_gg,
  filename = here("results/paper/figures/revision_cell_cycle_correction_LISI_scores.pdf"),
  width = 5,
  height = 11,
  units = "cm"
)



