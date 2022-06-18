# This script creates the raw data to create the main panels of figure 1


# Load packages
library(tidyverse)


# Read data
# rna_multiome_df <- readRDS("~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/final_clusters/dataframe_with_all_cells_20210929_with_multiome_myeloid.rds")
# cite_df <- readRDS("~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/final_clusters/figure_1_umap_df_cite_seq.rds")
# atac_df <- readRDS("~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/final_clusters/figure_1_umap_df_atac_seq.rds")
# donor_metadata <- read_csv("~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/data/tonsil_atlas_donor_metadata.csv")
rna_multiome_df <- readRDS(here::here("scRNA-seq/results/R_objects/final_clusters/dataframe_with_all_cells_20210929_with_multiome_myeloid.rds"))
cite_df <- readRDS(here::here("scRNA-seq/results/R_objects/final_clusters/figure_1_umap_df_cite_seq.rds"))
atac_df <- readRDS(here::here("scRNA-seq/results/R_objects/final_clusters/figure_1_umap_df_atac_seq.rds"))
donor_metadata <- read_csv(here::here("data/tonsil_atlas_donor_metadata.csv"))


# Clean, merge
rna_multiome_df <- left_join(
  rna_multiome_df,
  donor_metadata,
  by = c("donor_id", "hospital")
)
selected_cols <- c("barcode", "assay", "annotation_figure_1", "UMAP_1_level_1",
                   "UMAP_2_level_1", "donor_id", "sex", "age_group", "hospital")
figure_1_df <- rna_multiome_df %>%
  dplyr::select(all_of(selected_cols)) %>%
  dplyr::mutate(annotation_prob = NA)
figure_1_df$assay[figure_1_df$assay == "3P"] <- "scRNA-seq"
colnames(cite_df)[1:5] <- c("barcode", "annotation_figure_1", "annotation_prob",
                            "UMAP_1_level_1", "UMAP_2_level_1")
cite_df$assay <- "CITE-seq"
cite_df <- cite_df[, colnames(figure_1_df)]
colnames(atac_df)[1:5] <- c("barcode", "annotation_figure_1", "annotation_prob",
                            "UMAP_1_level_1", "UMAP_2_level_1")
atac_df$assay <- "scATAC-seq"
atac_df <- atac_df[, colnames(figure_1_df)]
figure_1_df <- bind_rows(list(figure_1_df, cite_df, atac_df))


# Save
# write_csv(
#   figure_1_df,
#   file = "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/tables/data_figures_paper/umaps_figure_1.csv"
# )
write_csv(
  figure_1_df,
  file = here::here("scRNA-seq/results/tables/data_figures_paper/umaps_figure_1.csv")
)