# This script creates the raw data to create the main panels of figure 1


# Load packages
library(tidyverse)
library(here)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Read data
rna_multiome_cite <- readRDS(path_to_save_cite_seq_df)
atac_multiome <- readRDS(path_to_save_atac_seq_df)


# UMAP figure 1
sel_cols_fig1 <- c("barcode", "UMAP_1_level_1", "UMAP_2_level_1",
                   "annotation_figure_1", "assay", "annotation_prob")
rna_multiome_cite_sub <- rna_multiome_cite[, sel_cols_fig1]
atac_multiome_sub <- atac_multiome[atac_multiome$type != "reference", sel_cols_fig1]
umap_fig_1_df <- bind_rows(rna_multiome_cite_sub, atac_multiome_sub)


# Save
umap_fig_1_df$assay <- case_when(
  umap_fig_1_df$assay == "3P" ~ "scRNA-seq",
  umap_fig_1_df$assay == "scATAC" ~ "scATAC-seq",
  umap_fig_1_df$assay == "multiome" ~ "Multiome",
  umap_fig_1_df$assay == "CITE-seq" ~ "CITE-seq"
)
rna_multiome_cite$assay <- case_when(
  rna_multiome_cite$assay == "3P" ~ "scRNA-seq",
  rna_multiome_cite$assay == "multiome" ~ "Multiome",
  rna_multiome_cite$assay == "CITE-seq" ~ "CITE-seq"
)
atac_multiome$assay <- case_when(
  atac_multiome$assay == "scATAC" ~ "scATAC-seq",
  atac_multiome$assay == "multiome" ~ "Multiome"
)
write_delim(umap_fig_1_df, path_to_df_umap_fig1, delim = ";", col_names = TRUE)
write_delim(rna_multiome_cite, path_to_qc_metrics_rna_multi_cite, delim = ";", col_names = TRUE)
write_delim(atac_multiome, path_to_qc_metrics_atac_multi, delim = ";", col_names = TRUE)