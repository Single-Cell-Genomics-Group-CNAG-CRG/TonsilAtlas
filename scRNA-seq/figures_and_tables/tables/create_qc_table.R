# This script creates a supplementary table with basic QC metrics per modality
# and donor


# Load package
library(tidyverse)
library(openxlsx)
library(here)
library(Seurat)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))


# Define parameters/variables
order_donors <- c("BCLL-2-T", "BCLL-6-T", "BCLL-8-T", "BCLL-9-T",
                  "BCLL-10-T", "BCLL-11-T", "BCLL-12-T", "BCLL-13-T",
                  "BCLL-14-T", "BCLL-15-T")


# Read data
rna_multiome_cite <- read_delim(
  path_to_qc_metrics_rna_multi_cite,
  delim = ";",
  col_names = TRUE
)
atac_multiome <- read_delim(
  path_to_qc_metrics_atac_multi,
  delim = ";",
  col_names = TRUE
)
spatial <- readRDS(path_to_spatial)


# Calculate basic QC RNA
qc_rna <- rna_multiome_cite %>%
  filter(assay == "scRNA-seq") %>%
  group_by(donor_id) %>%
  summarise(
    n_cells = n(),
    mean_n_counts = round(mean(nCount_RNA), 2),
    mean_n_features = round(mean(nFeature_RNA), 2),
    mean_pct_mitochondrial = round(mean(pct_mt), 2),
  ) %>%
  as.data.frame()
indx_rna <- match(order_donors[order_donors %in% qc_rna$donor_id], qc_rna$donor_id)
qc_rna <- qc_rna[indx_rna, ]


# Calculate basic QC CITE-seq
qc_cite <- rna_multiome_cite %>%
  filter(assay == "CITE-seq") %>%
  group_by(donor_id) %>%
  summarise(
    n_cells = n(),
    mean_n_counts_RNA = round(mean(nCount_RNA), 2),
    mean_n_features_RNA = round(mean(nFeature_RNA), 2),
    mean_n_counts_ADT = round(mean(nCount_ADT), 2),
    mean_n_features_ADT = round(mean(nFeature_ADT), 2),
    mean_pct_mitochondrial = round(mean(pct_mt), 2),
  ) %>%
  as.data.frame()
indx_cite <- match(order_donors[order_donors %in% qc_cite$donor_id], qc_cite$donor_id)
qc_cite <- qc_cite[indx_cite, ]


# Calculate basic QC scATAC-seq
qc_atac <- atac_multiome %>%
  filter(assay == "scATAC-seq") %>%
  group_by(donor_id) %>%
  summarise(
    n_cells = n(),
    mean_n_counts = round(mean(nFeature_peaks_macs), 2),
    mean_n_features = round(mean(nFeature_peaks_macs), 2),
    mean_TSS_enrichment = round(mean(TSS.enrichment), 2),
    mean_pct_reads_in_peaks = round(mean(pct_reads_in_peaks), 2),
  ) %>%
  as.data.frame()
indx_atac <- match(order_donors[order_donors %in% qc_atac$donor_id], qc_atac$donor_id)
qc_atac <- qc_atac[indx_atac, ]


# Calculate basic QC Multiome
qc_multiome1 <- rna_multiome_cite %>%
  filter(assay == "Multiome") %>%
  group_by(donor_id) %>%
  summarise(
    n_cells_RNA = n(),
    mean_n_counts_RNA = round(mean(nCount_RNA), 2),
    mean_n_features_RNA = round(mean(nFeature_RNA), 2),
    mean_pct_mitochondrial = round(mean(pct_mt), 2),
  ) %>%
  as.data.frame()
indx_multiome1 <- match(order_donors[order_donors %in% qc_multiome1$donor_id], qc_multiome1$donor_id)
qc_multiome1 <- qc_multiome1[indx_multiome1, ]

qc_multiome2 <- atac_multiome %>%
  filter(assay == "Multiome") %>%
  group_by(donor_id) %>%
  summarise(
    n_cells_ATAC = n(),
    mean_n_counts_ATAC =  round(mean(nFeature_peaks_macs), 2),
    mean_n_features_ATAC =  round(mean(nFeature_peaks_macs), 2),
    mean_TSS_enrichment_ATAC =  round(mean(TSS.enrichment), 2)
  ) %>%
  as.data.frame()
indx_multiome2 <- match(order_donors[order_donors %in% qc_multiome2$donor_id], qc_multiome2$donor_id)
qc_multiome2 <- qc_multiome2[indx_multiome2, ]

qc_multiome <- left_join(x = qc_multiome1, y = qc_multiome2, by = "donor_id")
qc_multiome <- qc_multiome[, c(1, 2, 6, 3, 4, 5, 7:ncol(qc_multiome))]


# Spatial
qc_spatial <- spatial@meta.data %>%
  group_by(donor_id) %>%
  summarise(
    n_spots = n(),
    mean_n_counts = round(mean(nCount_Spatial), 2),
    mean_n_features = round(mean(nFeature_Spatial), 2),
    mean_pct_mitochondrial = round(mean(pct_mt), 2),
  ) %>%
  as.data.frame()
indx_spatial <- match(order_donors[order_donors %in% qc_spatial$donor_id], qc_spatial$donor_id)
qc_spatial <- qc_spatial[indx_spatial, ]


# Save
qc_list <- list(
  "scRNA-seq" = qc_rna,
  "scATAC-seq" = qc_atac,
  "CITE-seq" = qc_cite,
  "Multiome" = qc_multiome,
  "Spatial" = qc_spatial
)
openxlsx::write.xlsx(x = qc_list, file = path_to_save_qc_table)
  
  