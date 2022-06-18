# This script plots the supplementary figure for cytotoxic cells

# Load packages
library(Seurat)
library(Signac)
library(scCustomize)
library(tidyverse)
library(ggrastr)
library(patchwork)
library(ggpubr)
library(here)
library(glue)


# Source utils
source(here("scRNA-seq/bin/utils_figure3.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))
path_to_wd <- here()


# Read data
cd8 <- readRDS(path_to_save_cd8)
ilc_nk <- readRDS(path_to_save_ilc_nk)
seurat_cite <- readRDS(path_to_save_cite_cytotoxic)
seurat_cite <- subset(seurat_cite, subset = annotation_prob >= 0.6)
atac <- readRDS(path_to_save_atac_cytotoxic)
DefaultAssay(seurat_cite) <- "ADT"


# Dot Plot (RNA) expression
goi_rna <- c("BACH2", "TCF7", "PLAC8", "EOMES", "ITGA1", "HLA-DRB1", "CD38",
             "ITGAE", "IFNG", "CD200", "CCL4", "CD99", "CX3CR1", "NKG7",
             "SELL", "IFNG-AS1", "CD160", "CD96", "GZMH", "IL4I1", "IL1R1",
             "AHR")
ilc_nk$UMAP_1_20220215 <- ilc_nk$UMAP_1_20220215 + 10
ilc_nk$UMAP_2_20220215 <- ilc_nk$UMAP_2_20220215 + 10
seurat_rna <- merge(x = cd8, y = ilc_nk)
seurat_rna$annotation_20220215[seurat_rna$annotation_20220215 == "CXCR6+ RM CD8 T"] <- "RM CD8 activated T"
seurat_rna$annotation_paper <- factor(
  seurat_rna$annotation_20220215,
  levels = rev(names(colors_rna))
)
seurat_rna <- ScaleData(seurat_rna, features = rownames(seurat_rna))
Idents(seurat_rna) <- "annotation_paper"
features <- unlist(goi_rna_supplement)
features <- features[!(features %in% goi_rna)]
features <- features[!duplicated(features)]
names(features) <- NULL
dot_plot_all <- DotPlot(seurat_rna, features = features) +
  # coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1)
  )

# CITE-seq
DefaultAssay(seurat_cite) <- "ADT"
seurat_cite$annotation_20220215[seurat_cite$annotation_20220215 == "CXCR6+ RM CD8 T"] <- "RM CD8 activated T"
seurat_cite$annotation_paper <- factor(
  seurat_cite$annotation_20220215,
  levels = names(colors_rna)
)
# goi_cite <- c("CD45RA", "CD45RO", "CD99", "CD28", "CD38", "CD279-(PD-1)",
#               "CD278-(ICOS)", "CD336-(NKp44)")
# goi_cite2 <- c("CD103-(Integrin-alphaE)", "CD54", "CD161",
#                "CD56-(NCAM)")
goi_cite <- c("CD95-(Fas)", "CD122-(IL-2Rbeta)", "CD183-(CXCR3)", "CD127-(IL-7Ralpha)", 
              "CD69", "CD244-(2B4)", "CD16")
cite_vlns <- map(goi_cite, function(x) {
  p <- VlnPlot(
    seurat_cite,
    x,
    group.by = "annotation_paper",
    cols = colors_rna,
    pt.size = 0
  ) &
    NoLegend() &
    stat_summary(
      fun = "median",
      geom = "point",
      color = "black",
      size = 0.75
    ) &
    ylab("Protein Expression") &
    theme(
      axis.title.x = element_blank(),
      axis.text = element_text(size = 6),
      axis.title.y = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5)
    )
  p
})
cite_vlns_arranged <- wrap_plots(cite_vlns, nrow = 3, ncol = 3)


# Save
ggsave(
  filename = here("results/paper/figures/supplementary_figure_cytotoxic.pdf"),
  plot = dot_plot_all,
  width = 21,
  height = 14,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_cytotoxic.pdf"),
  plot = cite_vlns_arranged,
  width = 20,
  height = 14,
  units = "cm"
)


