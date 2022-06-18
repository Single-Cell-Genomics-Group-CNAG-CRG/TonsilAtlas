# This script performs differential expression analysis between healthy
# and tumoral clusters.

# Load packages
library(Seurat)
library(ggrastr)
library(ggrepel)
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


# Differential expression analysis (Treg)
cd4_tonsil <- cd4_tonsil[, cd4_tonsil$assay == "multiome"] # to remove batch effects
cd4_tonsil$group <- "HD"
cd4_tonsil$cell_types <- cd4_tonsil$annotation_20220215
cd4_102$group <- "MCL"
cd4_102$cell_types <- cd4_102$annotation_SLOcatoR
cd4_413$group <- "MCL"
cd4_413$cell_types <- cd4_413$annotation_SLOcatoR
merged <- merge(x = cd4_tonsil, c(cd4_102, cd4_413))
Idents(merged) <- "cell_types"
treg_levels <- c("Eff-Tregs", "non-GC-Tf-regs", "GC-Tf-regs")
dea_tregs <- purrr::map(treg_levels, function(x) {
  seurat_obj <- merged[, merged$cell_types == x]
  Idents(seurat_obj) <- "group"
  dea <- FindMarkers(
    seurat_obj,
    ident.1 = "MCL",
    ident.2 = "HD",
    only.pos = FALSE,
    logfc.threshold = 0.3
  )
  dea$is_significant <- ifelse(
    dea$p_val_adj < 0.01 & abs(dea$pct.1 - dea$pct.2) > 0.25,
    "significant",
    "not significant"
  )
  dea$gene <- rownames(dea)
  dea
})
names(dea_tregs) <- treg_levels


# Volcano plot
top_genes <- c("MT2A", "MT1X", "MT1E", "STAT4", "TNFRSF9", "CXCR4", "IKZF2",
               "TSHZ2", "LAG3", "SLA", "IL7R", "MAF", "IL1R1")
dea_sub <- dea_tregs$`Eff-Tregs`[top_genes, ]
(volcano_gg <- dea_tregs$`Eff-Tregs` %>%
  ggplot(aes(avg_log2FC, -1 * log10(p_val_adj), color = is_significant)) +
    geom_point(size = 0.25) +
    geom_text_repel(data = dea_sub, aes(label = gene), color = "black",
                    size = 2) +
    scale_color_manual(
      values = c("gray", "forestgreen"),
      breaks = c("not significant", "significant"),
      labels = c("not sig.", "sig.")
    ) +
    labs(x =  "log2 (MCL / HD)", y = "-log10(adj. p-value)") +
    theme_classic()) +
    theme(
      legend.title = element_blank(),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6)) +
    guides(color = guide_legend(override.aes = list(size = 2)))



# Save
saveRDS(dea_tregs, here("data/raw_data_figures/dea_tregs_list.rds"))


