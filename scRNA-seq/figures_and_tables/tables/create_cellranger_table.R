# This script creates a table with the raw QC metrics coming from cellranger
# for every OMIC


# Load packages
library(tidyverse)
library(here)


# Define paths
path_to_rna_hash <- here("scRNA-seq/results/tables/cellranger_mapping/cellranger_mapping_metrics_hashed.csv")
path_to_rna_no_hash <- here("scRNA-seq/results/tables/cellranger_mapping/cellranger_mapping_metrics_not_hashed.csv")
path_to_atac <- here("scATAC-seq/results/tables/cellranger_mapping/cellranger_mapping_metrics_atac.csv")
path_to_multiome <- here("multiome/results/tables/cellranger_mapping/cellranger_mapping_metrics_multiome.csv")
path_to_cite <- here("CITE-seq/results/tables/cellranger_mapping/cellranger_mapping_metrics_cite_seq.csv")
path_to_spatial <- here("spatial_transcriptomics/02-QC/2020-09-22/R_objects_2020-09-22/Supplementary_QC_table_spatial.csv")
path_to_metadata_rna <- here("scRNA-seq/1-cellranger_mapping/data/tonsil_atlas_metadata.csv")
path_to_metadata_atac <- here("scATAC-seq/1-cellranger_mapping/data/tonsil_atlas_metadata_atac.csv")
path_to_metadata_multiome <- here("multiome/1-cellranger_mapping/data/tonsil_atlas_metadata_multiome.csv")
path_to_metadata_cite <- here("CITE-seq/1-cellranger_mapping/data/tonsil_atlas_metadata.csv")
path_to_save_cellranger_table <- here("results/paper/tables/supplementary_table_cellranger_metrics.xlsx")

# Read data
rna_hash <- read_csv(path_to_rna_hash)
rna_no_hash <- read_csv(path_to_rna_no_hash)
atac <- read_csv(path_to_atac)
multiome <- read_csv(path_to_multiome)
cite <- read_csv(path_to_cite)
# cite_rna <- read_csv(path_to_cite_rna)
# cite_adt <- read_csv(path_to_cite_adt)
spatial <- read_csv(path_to_spatial)
metadata_rna <- read_csv(path_to_metadata_rna)
metadata_atac <- read_csv(path_to_metadata_atac)
metadata_multiome <- read_csv(path_to_metadata_multiome)
metadata_cite <- read_csv(path_to_metadata_cite)


# RNA
metadata_rna <- metadata_rna %>%
  filter(type != "hashed_hto") %>%
  select("subproject", "gem_id", "donor_id")
rna <- bind_rows(rna_hash, rna_no_hash)
rna <- left_join(metadata_rna, rna, by = "gem_id")
rna <- rna %>%
  add_column(
    is_hashed = ifelse(!is.na(rna$`Antibody: Antibody Reads in Cells`), "yes", "no"),
    .after = "donor_id"
  )

# ATAC
metadata_atac <- metadata_atac[, c("subproject", "gem_id", "donor_id")]
atac <- left_join(metadata_atac, atac, by = "gem_id")


# Multiome
metadata_multiome <- metadata_multiome[metadata_multiome$type == "RNA", ]
metadata_multiome <- metadata_multiome[, c("subproject", "gem_id", "donor_id")]
colnames(multiome)[colnames(multiome) == "Sample ID"] <- "gem_id"
multiome <- left_join(metadata_multiome, multiome, by = "gem_id")


# CITE-seq
modalities <- c("Gene Expression", "Antibody Capture", "VDJ B", "VDJ T")
cite_dfs <- map(modalities, function(x) {
  df <- cite %>%
    filter(
      `Library or Sample` == "Sample",
      `Library Type` %in% x
    ) %>%
    pivot_wider(
      id_cols = c(subproject, gem_id),
      names_from = `Metric Name`,
      values_from = `Metric Value`
    )
  df
})
cols_exclude <- c("subproject", "Cells")
cite_dfs[[2]] <- cite_dfs[[2]][, !(colnames(cite_dfs[[2]]) %in% cols_exclude)]
colnames(cite_dfs[[1]])[colnames(cite_dfs[[1]]) == "Median UMI counts per cell"] <- "Median UMI counts per cell (RNA)"
colnames(cite_dfs[[2]])[colnames(cite_dfs[[2]]) == "Median UMI counts per cell"] <- "Median UMI counts per cell (ADT)"
cite_df <- left_join(cite_dfs[[1]], cite_dfs[[2]], by = "gem_id")
cite_df <- arrange(cite_df, subproject)
cite_dfs[[3]] <- cite_dfs[[3]][, colnames(cite_dfs[[3]]) != "subproject"]
cite_dfs[[4]] <- cite_dfs[[4]][, colnames(cite_dfs[[4]]) != "subproject"]
shared_cols <- intersect(colnames(cite_dfs[[3]]), colnames(cite_dfs[[4]]))
shared_cols <- shared_cols[shared_cols != "gem_id"]
colnames(cite_dfs[[3]])[colnames(cite_dfs[[3]]) %in% shared_cols] <- str_c(
  colnames(cite_dfs[[3]])[colnames(cite_dfs[[3]]) %in% shared_cols],
  "(BCR)",
  sep = " "
)
colnames(cite_dfs[[4]])[colnames(cite_dfs[[4]]) %in% shared_cols] <- str_c(
  colnames(cite_dfs[[4]])[colnames(cite_dfs[[4]]) %in% shared_cols],
  "(TCR)",
  sep = " "
)
vdj_df <- left_join(cite_dfs[[3]], cite_dfs[[4]], by = "gem_id")
cite_df_final <- left_join(cite_df, vdj_df, by = "gem_id")


# Save
qc_list <- list(
  "scRNA-seq" = rna,
  "scATAC-seq" = atac,
  "CITE-seq" = cite_df_final,
  "Multiome" = multiome,
  "Spatial" = spatial
)
openxlsx::write.xlsx(x = qc_list, file = path_to_save_cellranger_table)