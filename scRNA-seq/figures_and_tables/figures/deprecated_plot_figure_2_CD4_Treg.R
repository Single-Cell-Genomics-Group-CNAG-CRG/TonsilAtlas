# This script plots the main panels for Tregs in figure 2


# Load packages
library(Seurat)
library(Nebulosa)
library(tidyverse)
library(ggrastr)
library(here)


# Source utils
source(here("scRNA-seq/bin/utils_figure2.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
seurat_rna <- readRDS(path_to_save_cd4)
auc_mtx <-  read.csv(path_to_auc_mtx, row.names = 1, check.names = FALSE)
seurat_treg <- readRDS(here::here("scATAC-seq/results/R_objects/treg_atac.rds"))


# Violin plots activity
filt <- colnames(seurat_rna)[colnames(seurat_rna) %in% rownames(auc_mtx)]
auc_mtx <- auc_mtx[filt, c("FOXP3(+)", "RORC(+)"), drop = FALSE]
seurat_rna$FOXP3_activity <- NA
seurat_rna$FOXP3_activity[rownames(auc_mtx)] <- auc_mtx[, "FOXP3(+)"]
seurat_rna$RORC_activity <- NA
seurat_rna$RORC_activity[rownames(auc_mtx)] <- auc_mtx[, "RORC(+)"]
vars <- c("FOXP3_activity", "RORC_activity")
treg_levels <- c("Eff-Tregs", "non-GC-Tf-regs", "GC-Tf-regs")
vln_plots_activity <- map(vars, function(x) {
  p <- seurat_rna@meta.data %>%
    filter(annotation_20220215 %in% treg_levels & !is.na(seurat_rna@meta.data[[x]])) %>%
    ggplot(aes_string("annotation_20220215", x)) +
      geom_violin(fill = "red", alpha = 0.7) +
      labs(title = str_remove(x, "_activity"), y = "TF activity") +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 7.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      )
  p
})
# vln_plots_activity
# vln_plots_activity[[1]] <- vln_plots_activity[[1]] +
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
# vln_plots_activity[[1]] / vln_plots_activity[[2]]
vln_plots_activity <- vln_plots_activity &
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


# ATAC
motif_names <- unlist(seurat_treg@assays$peaks_level_5@motifs@motif.names)
seurat_treg$TCF7_motif <- seurat_treg@assays$chromvar@data["MA0769.2", ]
seurat_treg$LEF1_motif <- seurat_treg@assays$chromvar@data["MA0768.1", ]
seurat_treg$MAFG_motif <- seurat_treg@assays$chromvar@data["MA0659.2", ]
seurat_treg$RORC_motif <- seurat_treg@assays$chromvar@data["MA1151.1", ]
vars <- c("TCF7_motif", "LEF1_motif", "MAFG_motif", "RORC_motif")
vln_plots_access <- map(vars, function(x) {
  p <- seurat_treg@meta.data %>%
    ggplot(aes_string("annotation_20220215", x)) +
    geom_violin(fill = "green4", alpha = 0.7) +
    labs(title = str_remove(x, "_motif"), y = "Motif accessibility") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, size = 7.5)
    )
  p
})

vln_plots_access[[1]] <- vln_plots_access[[1]] +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
vln_plots_access[[2]] <- vln_plots_access[[2]] +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

(vln_plots_activity[[1]] / vln_plots_activity[[2]]) /
(vln_plots_access[[1]] / vln_plots_access[[2]] / vln_plots_access[[3]])
patchwork::wrap_plots(c(vln_plots_activity, vln_plots_access), ncol = 1)
