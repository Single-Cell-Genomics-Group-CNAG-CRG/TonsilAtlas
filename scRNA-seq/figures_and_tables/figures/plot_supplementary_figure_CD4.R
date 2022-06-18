# This script plots the supplementary ifgure associated with the CD4 T cells

# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(readxl)
library(ggrastr)
library(ggpubr)
library(scCustomize)
library(Nebulosa)
library(UCell)
library(here)


# Source utils
source(here("scRNA-seq/bin/utils_figure2.R"))
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
seurat <- readRDS(path_to_save_cd4)
seurat_atac <- readRDS(path_to_save_atac_cd4)
seurat_th <- readRDS(path_to_save_th)
seurat_cite <- readRDS(path_to_save_cite_cd4)
markers_pnas <- read_excel(
  here("data/raw_data_figures/pnas.1705551114.sd02.xlsx"),
  sheet = "cTfr vs eTreg",
  col_names = TRUE
)
auc_mtx <-  read.csv(path_to_auc_mtx, row.names = 1, check.names = FALSE)


# Dot plot Tcm
seurat$annotation_20220215 <- unfactor(seurat$annotation_20220215)
seurat$annotation_20220215[seurat$annotation_20220215 == "non-GC-Tf-regs"] <- "Eff-Tregs-IL32"
seurat$annotation_20220215[seurat$annotation_20220215 == "GC-Tf-regs"] <- "Tfr"
seurat$annotation_paper <- factor(
  seurat$annotation_20220215,
  levels = names(colors_rna)
)
goi_cm <- c("TCF7", "ICOS", "TIGIT", "CXCR5", "IL6ST", "PRDM1", "ANXA1", "S100A4", "ITGB1")
seurat_sub <- subset(seurat, idents = c("Naive", "CM Pre-non-Tfh", "CM PreTfh"))
dot_plot_n_cm <- DotPlot(seurat_sub, features = rev(goi_cm)) +
  coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1)
  )


# Accessibility score BCL6 (enhancer / promoter) across CD4 subtypes
seurat_atac$annotation_20220215 <- unfactor(seurat_atac$annotation_20220215)
seurat_atac$annotation_20220215[seurat_atac$annotation_20220215 == "non-GC-Tf-regs"] <- "Eff-Tregs-IL32"
seurat_atac$annotation_20220215[seurat_atac$annotation_20220215 == "GC-Tf-regs"] <- "Tfr"
seurat_atac$annotation_paper <- factor(
  seurat_atac$annotation_20220215,
  levels = names(colors_rna)
)
seurat_atac$annotation_paper <- factor(
  seurat_atac$annotation_20220215,
  levels = names(colors_rna)
)
box_acc <- purrr::map(c("BCL6_gene1", "BCL6_enhancer1"), function(x) {
  p <- seurat_atac@meta.data %>%
    ggplot(aes_string("annotation_paper", x, fill = "annotation_paper")) +
      geom_boxplot(outlier.size = 0.1) +
      scale_fill_manual(values = colors_rna, breaks = names(colors_rna)) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1,
                                   color = "black"),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 6)
      )
  p
})
box_acc[[1]] <- box_acc[[1]] +
  ylab("BCL6 ATAC Score (promoter)") +
  theme(axis.text.x = element_blank())
box_acc[[2]] <- box_acc[[2]] + ylab("BCL6 ATAC Score (enhancer)")
boxs <- box_acc[[1]] / box_acc[[2]]


# Violin plots key markers
goi_cd4_supp <- c("CCR7", "CXCR5", "PRDM1", "PDCD1", "ICOS", "TIGIT", "CD40LG", "IL21", "BTLA",
                  "IL12RB2", "TGFBR1", "CD200", "IL4", "CCR5", "IL7R",
                  "CD28", "ITGB1", "CTLA4", "FOXP3", "IL2RA")
# seurat$annotation_20220215 <- factor(
#   seurat$annotation_20220215,
#   levels = names(colors_rna)
# )
Idents(seurat) <- "annotation_paper"
dot_plot_all <- DotPlot(seurat, features = rev(goi_cd4_supp)) +
  coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1)
  )


# CITE-seq: CD28, CD29
seurat_cite$annotation_figure_2[seurat_cite$annotation_figure_2 == "non-GC-Tf-regs"] <- "Eff-Tregs-IL32"
seurat_cite$annotation_figure_2[seurat_cite$annotation_figure_2 == "GC-Tf-regs"] <- "Tfr"
seurat_cite$annotation_paper <- factor(
  seurat_cite$annotation_figure_2,
  levels = names(colors_rna)
)
DefaultAssay(seurat_cite) <- "ADT"
vlns_prot_l <- map(c("CD28", "CD29"), function(x) {
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
vlns_prot_l[[1]] <- vlns_prot_l[[1]] & theme(axis.text.x = element_blank())
vlns_prot <- vlns_prot_l[[1]] / vlns_prot_l[[2]]


# Nebulosa plots key markers
umaps_th_supp <- plot_density(
  seurat_th,
  goi_th_supp,
  reduction = "umap",
  size = 0.3,
  combine = FALSE
)
umaps_th_supp <- map(umaps_th_supp, function(p) {
  p + scale_color_viridis_c(option = "magma") &
    NoLegend() &
    theme(
      plot.title = element_text(size = 6, hjust = 0.5),
      axis.text = element_text(size = 0),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
})
umaps_th_supp <- ggarrange(plotlist = umaps_th_supp, ncol = 5, nrow = 2)


# Plot signature PNAS paper (Wing et al.)
markers_pnas <- markers_pnas %>%
  filter(p.value < 0.001 & abs(m.value) > 1) %>%
  filter(!(gene_id %in% c(cc.genes$s.genes, cc.genes$g2m.genes)))
signature_cTfr <- markers_pnas$gene_id[markers_pnas$m.value > 0]
signature_eTreg <- markers_pnas$gene_id[markers_pnas$m.value < 0]
treg_levels <- c("Eff-Tregs", "Eff-Tregs-IL32", "Tfr")
Idents(seurat) <- "annotation_paper"
seurat_treg <- subset(seurat, idents = treg_levels)
seurat_treg <- UCell::AddModuleScore_UCell(
  seurat_treg,
  features = list(eTreg_score = signature_eTreg, cTfr_score = signature_cTfr)
)
vlns_treg_l <- map(c("eTreg_score_UCell", "cTfr_score_UCell"), function(x) {
  p <- VlnPlot(
    seurat_treg,
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
    ylab("UCell Score") &
    theme(
      axis.title.x = element_blank(),
      axis.text = element_text(size = 6),
      axis.title.y = element_text(size = 6),
      plot.title = element_text(size = 8, hjust = 0.5)
    )
  p
})
vlns_treg_l[[1]] <- vlns_treg_l[[1]] + ggtitle("eTreg")
vlns_treg_l[[2]] <- vlns_treg_l[[2]] + ggtitle("cTfr")
vlns_treg_l[[1]] <- vlns_treg_l[[1]] & theme(axis.text.x = element_blank())
vlns_treg <- vlns_treg_l[[1]] / vlns_treg_l[[2]]


# Plot TF activity
filt <- colnames(seurat)[colnames(seurat) %in% rownames(auc_mtx)]
auc_mtx <- auc_mtx[filt, c("PRDM1(+)", "FOXP3(+)")]
seurat$FOXP3_activity <- NA
seurat$PRDM1_activity <- NA
seurat$FOXP3_activity[rownames(auc_mtx)] <- auc_mtx [, "FOXP3(+)"]
seurat$PRDM1_activity[rownames(auc_mtx)] <- auc_mtx[, "PRDM1(+)"]
vars1 <- c("FOXP3_activity", "PRDM1_activity")
seurat_rna_sub <- subset(
  seurat,
  cells = colnames(seurat)[!is.na(seurat$FOXP3_activity)]
)
vln_tfs_nm_cm <- seurat_rna_sub@meta.data %>%
  filter(annotation_paper %in% c("Naive", "CM Pre-non-Tfh", "CM PreTfh")) %>%
  ggplot(aes(annotation_paper, PRDM1_activity)) +
  geom_violin(fill = "red") +
  ylab("PRDM1 activity") +
  stat_summary(
    fun = "median",
    geom = "point",
    color = "black",
    size = 0.75
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, color = "black")
  )
vln_foxp3_treg <- seurat_rna_sub@meta.data %>%
  filter(annotation_paper %in% treg_levels) %>%
  ggplot(aes(annotation_paper, FOXP3_activity)) +
  geom_violin(fill = "red") +
  ylab("FOXP3 activity") +
  stat_summary(
    fun = "median",
    geom = "point",
    color = "black",
    size = 0.75
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, color = "black")
  )
vln_prdm1_treg <- seurat_rna_sub@meta.data %>%
  filter(annotation_paper %in% treg_levels) %>%
  ggplot(aes(annotation_paper, PRDM1_activity)) +
  geom_violin(fill = "red") +
  ylab("PRDM1 activity") +
  stat_summary(
    fun = "median",
    geom = "point",
    color = "black",
    size = 0.75
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, color = "black")
  )
vlns_tfs <- vln_foxp3_treg / vln_prdm1_treg

# Save
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_naive_CM.pdf"),
  dot_plot_n_cm,
  width = 10,
  height = 9,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_violins_BCL6_atac.pdf"),
  boxs,
  width = 10,
  height = 9,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_dotplot_all.pdf"),
  dot_plot_all,
  width = 14,
  height = 12,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_violins_protein.pdf"),
  vlns_prot,
  width = 7,
  height = 9.75,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_umaps_th_nebulosa.pdf"),
  umaps_th_supp,
  width = 9,
  height = 7,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_violins_treg.pdf"),
  vlns_treg,
  width = 4,
  height = 9.75,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_violins_PRDM1_CM_naive.pdf"),
  vln_tfs_nm_cm,
  width = 7,
  height = 4,
  units = "cm"
)
ggsave(
  filename = here("results/paper/figures/supplementary_figure_CD4_T_violins_PRDM1_FOXP3_Treg.pdf"),
  vlns_tfs,
  width = 7,
  height = 8,
  units = "cm"
)


