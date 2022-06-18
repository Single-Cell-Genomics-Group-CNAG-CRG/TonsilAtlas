# This script plots the main panels for Tregs
# Load packages
library(Seurat)
library(Signac)
library(scCustomize)
library(patchwork)
library(tidyverse)
library(here)
library(UCell)
library(openxlsx)
library(readxl)


# Souce utilities
source("scRNA-seq/bin/utils_final_clusters.R")
source("scRNA-seq/bin/utils_figure2.R")


# Read data
treg_levels <- c("Eff-Tregs", "non-GC-Tf-regs", "GC-Tf-regs")
seurat_rna <- readRDS(path_to_save_cd4)
seurat_rna <- subset(seurat_rna, idents = treg_levels)
path_to_atac_treg <- here("scATAC-seq/results/R_objects/Targetted_analysis/CD4_T/20220412_seurat_treg_atac.rds")
seurat_atac <- readRDS(path_to_atac_treg)
markers_pnas <- read_excel(
  here("data/raw_data_figures/pnas.1705551114.sd02.xlsx"),
  sheet = "cTfr vs eTreg",
  col_names = TRUE
)


# Dot plot RNA
goi_treg <- c("MAF", "LAG3", "CD2", "IL10", "IL1R1", "IL1R2", "IKZF1", "IKZF2", "IKZF3",
              "CADM1", "TNFRSF4", "TNFRSF9", "TNFRSF18",
              "LEF1", "TCF7", "RBMS3", "SESN3", "PDE3B")
dot_plot <- DotPlot(seurat_rna, features = rev(goi_treg)) +
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


# Motif plots
motif_treg <- c("MA0495.3", "MA1151.1", "MA0768.1", "MA0769.2")
names(motif_treg) <- c("MAFF", "RORC", "LEF1", "TCF7")
vln_plots_atac <- map(names(motif_treg), function(x) {
  p <- VlnPlot(
    seurat_atac,
    features = motif_treg[x],
    pt.size = 0,
    cols = colors_rna[treg_levels]
  ) +
    labs(title = x, y = "Accessibility") +
    scale_fill_manual(values = rep("chartreuse3", 3)) +
    stat_summary(fun = "mean", geom = "point", color = "black", size = 0.25) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 6),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 5.5),
      axis.text.y = element_text(size = 5)
    )
  p
})
vln_plots_atac[1:3] <- map(vln_plots_atac[1:3], function(p) {
  p + theme(axis.ticks.x = element_blank())
})
vln_plots_atac_arr <- wrap_plots(vln_plots_atac, ncol = 1)
motif_plots <- map(motif_treg, function(x) {
  p <- MotifPlot(seurat_atac, x, assay = "peaks_redefined")
  p <- p +
    theme(
      plot.title = element_text(size = 5),
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 5)
    )
  p
})
motif_plots_arr <- wrap_plots(motif_plots, ncol = 1)


# Save
# fig <- dot_plot | motif_plots_arr | vln_plots_atac
ggsave(
  here("results/paper/figures/figure_2_treg_dotplot.pdf"),
  dot_plot,
  height = 8.5,
  width = 5,
  units = "cm"
)
ggsave(
  here("results/paper/figures/figure_2_treg_dotplot_legend.pdf"),
  legend_dot_plot,
  height = 8.5,
  width = 5,
  units = "cm"
)
ggsave(
  here("results/paper/figures/figure_2_treg_motifs_logos.pdf"),
  motif_plots_arr,
  height = 10.5,
  width = 4.5,
  units = "cm"
)
ggsave(
  here("results/paper/figures/figure_2_treg_motifs_violin.pdf"),
  vln_plots_atac_arr,
  height = 9.5,
  width = 4,
  units = "cm"
)
