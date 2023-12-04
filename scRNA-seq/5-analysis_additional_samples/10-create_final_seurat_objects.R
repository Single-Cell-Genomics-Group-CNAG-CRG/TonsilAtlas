# This script reads all the annotated Seurat objects and fetches the metadata.
# Subsequently, it merges all data frames and subsets cells from the main
# Seurat object to discard doublets/poor-quality cells and adds the annotation


# Load packages
print("Loading packages...")
library(Seurat)
library(harmony)
library(tidyverse)
library(here)
library(glue)


# Load functions
print("Loading functions...")
source(here("scRNA-seq/bin/SLOcatoR_functions.R"))
source(here("scRNA-seq/bin/utils.R"))


# Read data
print("Reading data...")
donor_metadata <- read_csv(here("data/tonsil_atlas_donor_metadata.csv"))
directory <- here("scRNA-seq/results/R_objects/seurat_objects_revision")
files <- c(
  myeloid = "5.1-myeloid_seurat_object_discovery_validation_cohorts_annotation.rds",
  epithelial = "5.2-epithelial_seurat_object_discovery_validation_cohorts_annotation.rds",
  PDC = "5.3-PDC_seurat_object_discovery_validation_cohorts_annotation.rds",
  FDC = "5.4-FDC_seurat_object_discovery_validation_cohorts_annotation.rds",
  preB = "4.5-preBC_seurat_object_discovery_validation_cohorts.rds",
  preT = "4.6-preTC_seurat_object_discovery_validation_cohorts.rds",
  Cytotoxic = "5.7-Cytotoxic_seurat_object_discovery_validation_cohorts_annotation.rds",
  CD4_T = "5.8-CD4_seurat_object_discovery_validation_cohorts_annotation.rds",
  PC = "5.9-pc_seurat_object_discovery_validation_cohorts_annotation.rds",
  NBC_MBC = "5.10-NBC_MBC_seurat_object_discovery_validation_cohorts_annotation.rds",
  GCBC = "5.11-GCBC_seurat_object_discovery_validation_cohorts_annotation.rds"
)
cell_types <- names(files)
files <- str_c(directory, files, sep = "/")
seurat_list <- map(files, readRDS)
names(seurat_list) <- cell_types


# Update metadata
excluded_cols <- c("hospital", "sex", "age", "age_group", "cause_for_tonsillectomy", "cohort_type", "comments")
selected_cols <- c(
  "barcode", "donor_id", "gem_id", "library_name", "assay", "sex", "age",
  "age_group", "hospital", "cohort_type", "cause_for_tonsillectomy", "is_hashed", "preservation",
  "nCount_RNA", "nFeature_RNA", "pct_mt", "pct_ribosomal", "pDNN_hashing", "pDNN_scrublet", "pDNN_union", "scrublet_doublet_scores",
  "S.Score", "G2M.Score", "Phase", "scrublet_predicted_doublet", "doublet_score_scDblFinder",
  "annotation_level_1", "annotation_level_1_probability", "annotation_figure_1",
  "annotation_20220215", "annotation_20220619", "annotation_20230508",
  "annotation_20230508_probability", "UMAP_1_level_1", "UMAP_2_level_1",
  "UMAP_1_20220215", "UMAP_2_20220215", "UMAP_1_20230508", "UMAP_2_20230508"
)
seurat_list$preB$annotation_20220619 <- "preB"
seurat_list$preB$annotation_20230508 <- "preB"
seurat_list$preB$annotation_20230508_probability <- NA
seurat_list$preB$UMAP_1_20230508 <- NA
seurat_list$preB$UMAP_2_20230508 <- NA

seurat_list$preT$annotation_20220619 <- "preT"
seurat_list$preT$annotation_20230508 <- "preT"
seurat_list$preT$annotation_20230508_probability <- NA
seurat_list$preT$UMAP_1_20230508 <- NA
seurat_list$preT$UMAP_2_20230508 <- NA

seurat_list <- map(seurat_list, \(seurat_obj) {
  old_metadata <- seurat_obj@meta.data
  old_metadata <- old_metadata[, which(!colnames(old_metadata) %in% excluded_cols)]
  new_metadata <- as.data.frame(left_join(old_metadata, donor_metadata, by = "donor_id"))
  rownames(new_metadata) <- new_metadata$barcode
  new_metadata$UMAP_1_20230508 <- Embeddings(seurat_obj, "umap")[, "UMAP_1"]
  new_metadata$UMAP_2_20230508 <- Embeddings(seurat_obj, "umap")[, "UMAP_2"]
  new_metadata <- new_metadata[, selected_cols]
  seurat_obj@meta.data <- new_metadata
  seurat_obj
})


# Polish annotation
seurat_list$myeloid$annotation_20230508 <- as.character(seurat_list$myeloid$annotation_20230508)
seurat_list$myeloid$annotation_20230508[seurat_list$myeloid$annotation_20230508 == "MMP Slancytes"] <- "MMP Slan-like"
seurat_list$myeloid$annotation_20230508[seurat_list$myeloid$annotation_20230508 == "C1Q Slancytes"] <- "C1Q Slan-like"
seurat_list$myeloid$annotation_20230508[seurat_list$myeloid$annotation_20230508 == "SELENOP Slancytes"] <- "SELENOP Slan-like"
seurat_list$myeloid$annotation_20230508[seurat_list$myeloid$annotation_20230508 == "ITGAX Slancytes"] <- "ITGAX Slan-like"
new_myeloid_levels <- c("DC1 precursor", "DC1 mature", "DC2", "DC3", "DC4",
                        "DC5", "IL7R DC", "aDC1", "aDC2", "aDC3", "M1 Macrophages",
                        "Monocytes", "Mast", "Neutrophils", "Cycling",
                        "MMP Slan-like", "C1Q Slan-like", "SELENOP Slan-like",
                        "ITGAX Slan-like")
seurat_list$myeloid$annotation_20230508 <- factor(
  seurat_list$myeloid$annotation_20230508,
  new_myeloid_levels
)
seurat_list$Cytotoxic$annotation_20230508 <- as.character(seurat_list$Cytotoxic$annotation_20230508)
seurat_list$Cytotoxic$annotation_20230508[seurat_list$Cytotoxic$annotation_20230508 == "CD16-CD56- NK"] <- "CD16-CD56dim NK"
new_cytotoxic_levels <- c("Naive CD8 T", "SCM CD8 T", "CM CD8 T", "RM CD8 T",
                          "RM CD8 activated T", "CD8 Tf", "DC recruiters CD8 T",
                          "IFN+ CD8 T", "EM CD8 T", "ZNF683+ CD8 T", "TCRVδ+ gd T",
                          "MAIT/CD161+TRDV2+ gd T-cells", "CD16-CD56+ NK", "CD16-CD56dim NK",
                          "CD16+CD56- NK", "ILC1", "NKp44+ ILC3", "NKp44- ILC3",
                          "DN")
seurat_list$Cytotoxic$annotation_20230508 <- factor(
  seurat_list$Cytotoxic$annotation_20230508,
  new_cytotoxic_levels
)
Idents(seurat_list$Cytotoxic) <- "annotation_20230508"
Idents(seurat_list$myeloid) <- "annotation_20230508"
Idents(seurat_list$PDC) <- "annotation_20230508"

# Save Seurat objects
path_save_seurat <- here("scRNA-seq/results/R_objects/final_clusters/20230911")
for (x in names(seurat_list)) {
  print(x)
  saveRDS(seurat_list[[x]], glue("{path_save_seurat}/20230911_{x}_seurat_obj.rds"))
}


# Get dataframes
print("Getting metadata dataframes for each Seurat object...")
metadata_dfs <- map(seurat_list, \(seurat_obj) {
  df <- seurat_obj@meta.data
  df
})
metadata_df <- bind_rows(metadata_dfs)


# Include unannotated multiome cells from underrepresented cell types
print("Including unannotated multiome cells from underrepresented cell types...")
epithelial <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/4.2-epithelial_seurat_object_discovery_validation_cohorts.rds"))
epithelial <- updateAnnotation(epithelial)
epithelial_multiome <- epithelial@meta.data[epithelial$assay == "multiome", ]
fdc <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/4.4-FDC_seurat_object_discovery_validation_cohorts.rds"))
fdc <- updateAnnotation(fdc)
fdc_multiome <- fdc@meta.data[fdc$assay == "multiome", ]
pdc <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/4.3-PDC_seurat_object_discovery_validation_cohorts.rds"))
pdc <- updateAnnotation(pdc)
pdc_multiome <- pdc@meta.data[pdc$assay == "multiome", ]
multiome_df <- bind_rows(epithelial_multiome, fdc_multiome, pdc_multiome)
multiome_df$annotation_20230508 <- "unannotated"
multiome_df$annotation_20230508_probability <- NA
multiome_df$UMAP_1_20230508 <- NA
multiome_df$UMAP_2_20230508 <- NA
multiome_df <- multiome_df[, selected_cols]

# Merge, subset and add metadata
print("Arranging metadata...")
if (all(colnames(metadata_df) == colnames(multiome_df))) {
  metadata_df <- bind_rows(metadata_df, multiome_df)
} 
seurat <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/3.1-seurat_object_discovery_validation_cohorts_integrated.rds"))
seurat <- seurat[, metadata_df$barcode]
if (all(colnames(seurat) == metadata_df$barcode)) {
  seurat@meta.data <- metadata_df
}


# Reprocess
print("Processing whole object...")
seurat$type <- str_c(seurat$cohort_type, seurat$assay, sep = "_")
shared_hvg <- find_assay_specific_features(seurat, assay_var = "type")
print(glue("The number of shared HVG is: {length(shared_hvg)}"))
seurat <- integrate_assays(
  seurat,
  assay_var = "type",
  shared_hvg = shared_hvg
)
seurat <- RunUMAP(seurat, dims = 1:30, reduction = "harmony")
seurat$UMAP_1_20230508 <- Embeddings(seurat, "umap")[, "UMAP_1"]
seurat$UMAP_2_20230508 <- Embeddings(seurat, "umap")[, "UMAP_2"]
new_levels <- c("NBC_MBC", "GCBC", "PC", "CD4_T", "Cytotoxic",
                "myeloid", "FDC", "PDC", "epithelial", "preBC", "preTC")
seurat$annotation_level_1 <- factor(seurat$annotation_level_1, new_levels)
Idents(seurat) <- "annotation_level_1"

# Save
print("Saving...")
path_save <- here("scRNA-seq/results/R_objects/final_clusters/20230911/20230911_tonsil_atlas_rna_seurat_obj.rds")
#path_save <- here("scRNA-seq/results/R_objects/seurat_objects_revision/7-20230529_seurat_rna_discovery_validation_cohorts.rds")
saveRDS(seurat, path_save)
