# This script prepares the 10X scATAC-seq object to share in Zenodo

# For ATAC, we need to:
#     1. Select all cells that we will include. That's a mix of the all cells from final scATAC-seq object and multiome cells from final multiome object.
#     2. Rename cell barcodes for scATAC-seq: include prefix
#     3. Identify fragment files, only keep fragments from cells in the object. Save paths to fragment files to include in object.
#     4. Call peaks per cell type. The purpose is to identify peaks that are specific to all cell types, even rare ones. We might need to merge cell types to ensure we dont have too few cells
#     5. Create accessibility matrix based on the previous peaks (4) and the fragments file (3)
#     6. Normalize and perform dimensionality reduction


# Load packages
library(Seurat)
library(Signac)
library(caret)
library(class)
library(tidyverse)
library(here)
library(glue)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(harmony)
plan("multicore", workers = 20)
options(future.globals.maxSize = 200000 * 1024^2) # for 100 Gb RAM


# Load object that contains RNA+multiome, subset to keep multiome only
seurat_atac <- readRDS(here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_tonsil_atlas_atac_seurat_obj.rds"))
seurat_multiome <- readRDS(here("scRNA-seq/results/R_objects/final_clusters/20230911/20230911_tonsil_atlas_multiome_seurat_obj.rds"))
donor_metadata <- read_csv(here("data/tonsil_atlas_donor_metadata.csv"))

  
# Select cells
seurat_atac <- seurat_atac[, seurat_atac$assay == "scATAC"]
new_names_atac <- str_c(
  seurat_atac$gem_id,
  "_",
  str_sub(colnames(seurat_atac), 1, 17),
  "1",
  sep = ""
)
seurat_atac <- RenameCells(seurat_atac, new.names = new_names_atac)
seurat_atac$barcode <- colnames(seurat_atac)
selected_cells <- list()
for (gem_id in unique(seurat_atac$gem_id)) {
  selected_cells[[gem_id]] <- colnames(seurat_atac)[seurat_atac$gem_id == gem_id]
}
for (gem_id in unique(seurat_multiome$gem_id)) {
  selected_cells[[gem_id]] <- colnames(seurat_multiome)[seurat_multiome$gem_id == gem_id]
}


# Set fragments files
fragment_objs <- map(names(selected_cells), \(gem_id) {
  path_frags <- glue("./fragments_atac/{gem_id}_atac_fragments_with_prefix.tsv.gz")
  cells <- selected_cells[[gem_id]]
  names(cells) <- cells
  frags_obj <- CreateFragmentObject(path = path_frags, cells = cells)
  frags_obj
})
names(fragment_objs) <- names(selected_cells)


# Deprecatng this code block, instead we use the same peaks from multiome, which are already cell type specific on bonda-fide clusters
# Call peaks 
# peaks_list <- map(fragment_objs, \(frags_obj) {
#   peaks <- CallPeaks(
#     object = frags_obj,
#     macs2.path = "~/anaconda3/envs/richter/bin/macs2"
#   )
# })
# combined_peaks <- GenomicRanges::reduce(unlist(GRangesList(peaks_list)))
# #saveRDS(combined_peaks, here("scRNA-seq/results/R_objects/combined_peaks_scATAC_zenodo.rds"))

# Create accessibility matrix
peaks <- granges(seurat_multiome[["ATAC"]])
acc_mat <- FeatureMatrix(
  features = peaks,
  fragment_objs,
  cells = unlist(selected_cells)
)


# Create Chromatin Assay
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
chrom_assay <- CreateChromatinAssay(
  counts = acc_mat,
  sep = c("-", "-"),
  fragments = fragment_objs,
  annotation = annotation
)
# saveRDS(chrom_assay, here("scRNA-seq/results/R_objects/chromatin_assay_scATAC_zenodo.rds"))


# Arrange metadata
metadata_multiome <- seurat_multiome@meta.data
# rm(seurat_multiome)
metadata_multiome$is_facs_sorted <- "NA"
excluded_cols <- c("UMAP_1_level_1", "UMAP_2_level_1", "UMAP_1_20220215",
                   "UMAP_2_20220215", "UMAP_1_20230508", "UMAP_2_20230508", "type",                           
                   "nucleosome_signal", "nucleosome_percentile", "TSS.enrichment",
                   "TSS.percentile", "RNA.weight", "ATAC.weight")
metadata_multiome <- metadata_multiome[, !(colnames(metadata_multiome) %in% excluded_cols)]
metadata_multiome$scrublet_predicted_doublet <- as.logical(metadata_multiome$scrublet_predicted_doublet)

metadata_atac <- seurat_atac@meta.data[colnames(seurat_atac) %in% colnames(chrom_assay), ]
selected_cols <- c("barcode", "donor_id", "gem_id", "library_name", "assay", "is_facs_sorted",
                   "sex", "age", "age_group", "hospital", "nCount_ATAC", "nFeature_ATAC",
                   "annotation_level_1", "annotation_figure_1", "scrublet_doublet_scores_atac",
                   "scrublet_predicted_doublet_atac")
metadata_atac <- metadata_atac[, selected_cols]
metadata_atac$cohort_type <- "NA"
metadata_atac$preservation <- "frozen"
metadata_atac$is_hashed <- "not_hashed"
metadata_atac$nCount_RNA <- NA
metadata_atac$nFeature_RNA <- NA
metadata_atac$pct_mt <- NA
colnames(metadata_atac)[colnames(metadata_atac) == "scrublet_doublet_scores_atac"] <- "scrublet_doublet_scores"
colnames(metadata_atac)[colnames(metadata_atac) == "scrublet_predicted_doublet_atac"] <- "scrublet_predicted_doublet"
metadata_atac$pct_ribosomal <- NA
metadata_atac$pDNN_hashing <- NA
metadata_atac$pDNN_union <- NA
metadata_atac$pDNN_scrublet <- NA
metadata_atac$S.Score <- NA
metadata_atac$G2M.Score <- NA
metadata_atac$Phase <- NA
metadata_atac$doublet_score_scDblFinder <- NA
joined_metadata <- as.data.frame(left_join(metadata_atac, donor_metadata, by = "donor_id"))
rownames(joined_metadata) <- joined_metadata$barcode
if (all(rownames(joined_metadata) == rownames(metadata_atac))){
  metadata_atac$cause_for_tonsillectomy <- joined_metadata$cause_for_tonsillectomy
}
metadata_atac$age_group[metadata_atac$age_group == "kid"] <- "child" # TODO: change "kid" for "child" in multiome, calculate QC metrics in multiome, handle annotations scATAC-seq data
metadata_atac$annotation_level_1_probability <- NA
metadata_atac$annotation_20220215 <- NA
metadata_atac$annotation_20220619 <- NA
metadata_atac$annotation_20230508 <- NA
metadata_atac$annotation_20230508_probability <- NA


# Create Seurat object
selected_cols <- c("barcode", "donor_id", "gem_id", "library_name", "assay",
                   "is_facs_sorted", "is_hashed", "sex", "age", "age_group", "hospital",
                   "cause_for_tonsillectomy", "nCount_ATAC", "nFeature_ATAC",
                   "annotation_level_1", "annotation_figure_1", "scrublet_doublet_scores",
                   "scrublet_predicted_doublet", "cohort_type", "preservation",
                   "pct_mt", "pct_ribosomal", "pDNN_hashing", "pDNN_union", "pDNN_scrublet",
                   "S.Score", "G2M.Score", "Phase", "doublet_score_scDblFinder", "annotation_level_1_probability",
                   "annotation_20220215", "annotation_20220619", "annotation_20230508", "annotation_20230508_probability")
if (all(selected_cols %in% colnames(metadata_multiome)) & all(selected_cols %in% colnames(metadata_atac))) {
  metadata_atac <- metadata_atac[, selected_cols]
  metadata_multiome <- metadata_multiome[, selected_cols]
  metadata_all <- as.data.frame(bind_rows(metadata_atac, metadata_multiome))
} else {
  stop("Selected columns are not in the dataframe!")
}
rownames(metadata_all) <- metadata_all$barcode
metadata_all <- metadata_all[colnames(chrom_assay), ]
seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = metadata_all
)


# Calculate QC metrics ATAC
# Calculate ATAC QC variables
seurat <- NucleosomeSignal(object = seurat)
seurat <- TSSEnrichment(object = seurat, fast = FALSE)


# Process ATAC-seq
seurat <- seurat %>%
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = "q0") %>%
  RunSVD()
seurat <- RunHarmony(
  seurat,
  reduction = "lsi",
  assay.use = "ATAC",
  reduction.save = "harmony",
  project.dim = FALSE,
  group.by.vars = "donor_id"
)
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 2:40)
#saveRDS(seurat, here("scRNA-seq/results/R_objects/seurat_atac_integrated_zenodo_v2.rds"))


# SLOcatoR
#source(here("scRNA-seq/bin/SLOcatoR_functions.R"))
#data_sets <- split_training_and_test_sets(
 # seurat,
 # split_var = "assay",
 # referece_label = "multiome",
  #query_label = "scATAC",
  #reduction = "harmony",
  #n_dims = 40
#)
#annotation_query_df <- transfer_label(
#  seurat_obj = seurat,
#  training_set = data_sets$training_set,
#  test_set = data_sets$test_set,
#  k = 8,
#  response_var = "annotation_20230508"
#)
#seurat$annotation_20230508[annotation_query_df$query_cells] <- annotation_query_df$annotation
#seurat$annotation_20230508_probability[annotation_query_df$query_cells] <- annotation_query_df$annotation_prob
# saveRDS(seurat, here("scRNA-seq/results/R_objects/seurat_atac_integrated_zenodo_v2_annotated.rds"))


# Annotations
#rm(seurat_atac, seurat_multiome)
#paths_objects_2022 <- c(
#  "NBC_MBC" = here("scATAC-seq/results/R_objects/Targetted_analysis/NBC_MBC/NBC_MBC_integrated_level_3.rds"),
#  "GCBC" = here("scATAC-seq/results/R_objects/Targetted_analysis/GCBC/GCBC_chromVar_CISBP_level_4.rds"),
#  "PC" = here("scATAC-seq/results/R_objects/Targetted_analysis/PC/05.PC_chromVar_CISBP_level_5.rds"),
#  "CD4_T" = here("scATAC-seq/results/R_objects/scATAC_CD4_T/files_plots/CD4T_integration_peak_calling_level_5.rds"),
#  "cytotoxic" = here("scATAC-seq/results/R_objects/Targetted_analysis/Cytotoxic/Cytotoxic_integrated_level_3.rds")
#)
#seurat_l_atac <- map(paths_objects_2022, readRDS)


# Save
saveRDS(seurat, here("scRNA-seq/results/R_objects/final_clusters/20230911/20230911_tonsil_atlas_atac_seurat_obj.rds"))
