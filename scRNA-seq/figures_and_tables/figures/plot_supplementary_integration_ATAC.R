# This script plots the supplementary figure containing:
# 1. UMAP scATAC-seq/multiome before and after integrating with harmony
# 2. LISI scores scATAC-seq/multiome before and after integrating with harmony


# Load packages
library(Seurat)
library(Signac)
library(lisi)
library(harmony)
library(patchwork)
library(ggrastr)
library(ggpubr)
library(here)
library(tidyverse)



# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure1.R"))
path_to_lisi_scores <- here("scATAC-seq/2-QC/5-batch_effect_correction/2-data_integration_multiome/tmp/lisi_scores.rds")


# Define parameters
seurat_atac <- readRDS(path_to_save_tonsil_atac)
lisi_scores <- readRDS(path_to_lisi_scores)


# UMAPs
cols <- cols_assays_df$color[cols_assays_df$assay %in% c("scATAC-seq", "Multiome")]
names(cols) <- cols_assays_df$assay[cols_assays_df$assay %in% c("scATAC-seq", "Multiome")]
umap_df_integr <- as.data.frame(Embeddings(seurat_atac, "umap"))
umap_df_integr$assay <- seurat_atac$assay
umap_df_integr$barcode <- rownames(umap_df_integr)
seurat_atac <- RunUMAP(
  seurat_atac,
  dims = 2:40,
  reduction = "lsi",
  reduction.name = "lsi_UMAP"
)
umap_df_not_integr <- as.data.frame(Embeddings(seurat_atac, "lsi_UMAP"))
colnames(umap_df_not_integr) <- c("UMAP_1", "UMAP_2")
umap_df_not_integr$assay <- seurat_atac$assay
umap_df_not_integr$barcode <- rownames(umap_df_not_integr)
umap_df <- bind_rows(
  list(unintegrated = umap_df_not_integr, integrated = umap_df_integr),
  .id = "is_integrated"
)
umap_df$assay <- case_when(
  umap_df$assay == "scATAC" ~ "scATAC-seq",
  umap_df$assay == "multiome" ~ "Multiome"
)
umap_df$assay <- factor(umap_df$assay, levels = names(cols))
umap_df$is_integrated <- factor(
  umap_df$is_integrated,
  levels = c("unintegrated", "integrated")
)
umaps_batch <- umap_df %>%
  ggplot(aes(UMAP_1, UMAP_2, color = assay)) +
  geom_point(shape = ".", alpha = 0.6) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  facet_grid(assay ~ is_integrated) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 6),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(size = 0.25)
  )
# umaps_batch <- rasterize(umaps_batch, dpi = 300)


# Plot lisi scores
lisi_scores$barcode <- rownames(lisi_scores)
lisi_scores$barcode <- str_remove(lisi_scores$barcode, "\\.\\.\\..*$")
lisi_scores <- lisi_scores[lisi_scores$barcode %in% colnames(seurat_atac), ]
lisi_vars <- c("lisi_sex", "lisi_age_group", "lisi_hospital", "lisi_assay")
confounders_vars <- c("sex", "age_group", "hospital", "assay")
confounders_df <- seurat_atac@meta.data[, c(confounders_vars, "barcode")]
colnames(lisi_scores) <- str_c("lisi", colnames(lisi_scores), sep = "_")
colnames(lisi_scores)[colnames(lisi_scores) == "lisi_is_integrated"] <- "is_integrated"
colnames(lisi_scores)[colnames(lisi_scores) == "lisi_barcode"] <- "barcode"
lisi_scores <- left_join(lisi_scores, confounders_df, by = "barcode")
cols_integration <- c(unintegrated = "gray60", integrated = "limegreen")
lisi_scores$is_integrated <- factor(
  lisi_scores$is_integrated,
  levels = names(cols_integration)
)
lisi_scores$age_group <- factor(
  lisi_scores$age_group,
  levels = c("kid", "young_adult", "old_adult")
)
levels(lisi_scores$age_group) <- c("kid", "young adult", "old adult")
lisi_scores$hospital <- factor(lisi_scores$hospital, levels = c("CIMA", "Clinic"))
levels(lisi_scores$hospital) <- c("Pamplona", "Barcelona")
lisi_scores$assay <- case_when(
  lisi_scores$assay == "scATAC" ~ "scATAC-seq",
  lisi_scores$assay == "multiome" ~ "Multiome"
)
lisi_scores$assay <- factor(
  lisi_scores$assay,
  levels = c("scATAC-seq", "Multiome")
)
boxplots_lisi <- purrr::map2(confounders_vars, lisi_vars, function(x, y) {
  p <- ggplot(lisi_scores, aes_string(x, y, fill = "is_integrated")) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(
      values = cols_integration,
      breaks = names(cols_integration)
    ) +
    labs(title = x, y = "LISI") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 7),
      axis.text = element_text(size = 6, color = "black"),
      axis.title.y = element_text(size = 6), 
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      axis.line = element_line(size = 0.25)
    )
  p
})
boxplots_arranged <- ggarrange(
  plotlist = boxplots_lisi,
  nrow = 2,
  ncol = 2,
  common.legend = TRUE
)


# Arrange figure
fig <- umaps_batch | boxplots_arranged


# Save
ggsave(
  plot = fig,
  filename = here("results/paper/figures/supplementary_figure_integration_ATAC.png"),
  width = 20,
  height = 10.5,
  dpi = 300,
  units = "cm"
)

