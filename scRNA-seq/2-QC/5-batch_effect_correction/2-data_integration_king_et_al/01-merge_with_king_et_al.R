# This script merges the Seurat object that does not contain doublets with the
# one obtain from Hamish King et al. (doi: 10.1126/sciimmunol.abe6291)


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_data <- here::here("scRNA-seq/results/R_objects/seurat_without_doublets_missing_reintegration.rds")
path_to_king_et_al <- here::here("data/King_et_al/CompleteIntegrated_scRNA_SeuratObject.rds")
path_to_metadata <- here::here("data/King_et_al/E-MTAB-8999.sdrf.txt")
path_to_save <- here::here("scRNA-seq/results/R_objects/seurat_merged_with_king_et_al.rds")


# Read data
tonsil_cnag <- readRDS(path_to_data)
tonsil_king <- readRDS(path_to_king_et_al)
metadata <- read_tsv(path_to_metadata, col_names = TRUE)


# Eliminate assays
DefaultAssay(tonsil_king) <- "RNA"
tonsil_king[["integrated"]] <- NULL


# Eliminate reductions
tonsil_cnag@reductions$pca <- NULL
tonsil_cnag@reductions$harmony <- NULL
tonsil_cnag@reductions$umap <- NULL


# Subset to only keep full tonsils (not enriched)
selected_samples <- str_subset(unique(tonsil_king$Sample), "_Total$")
all_cells <- colnames(tonsil_king)
selected_cells <- all_cells[tonsil_king$Sample %in% selected_samples]
tonsil_king <- subset(tonsil_king, cells = selected_cells)


# Homogenize QC metrics naming
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "percent.mito"] <- "pct_mt"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "percent.ribo"] <- "pct_ribosomal"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "ScrubletPrediction"] <- "scrublet_predicted_doublet"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "DoubletFinderPrediction"] <- "doublet_finder_predicted_doublet"
tonsil_cnag$doublet_finder_predicted_doublet <- NA
tonsil_king$HTO_classification.global <- NA
tonsil_king$pDNN_hashing <- NA
tonsil_king$pDNN_scrublet <- NA
tonsil_king$pDNN_union <- NA
tonsil_king$scrublet_doublet_scores <- NA
tonsil_king$scrublet_doublet_scores_scaled <- NA
tonsil_king$has_high_lib_size <- NA
tonsil_cnag$prePB_signature1 <- NULL
tonsil_cnag$apoptosis_signature2 <- NULL


# Define age, age group, sex and sampling center
metadata <- metadata %>% 
  filter(`Source Name` %in% selected_samples) %>%
  select("Source Name", "Characteristics[individual]", "Characteristics[sex]",
         "Characteristics[age]") %>%
  group_by(`Source Name`) %>%
  filter(row_number(`Source Name`) == 1) %>%
  magrittr::set_colnames(c("sample", "donor_id", "sex", "age"))
sex_vec <- metadata$sex
age_vec <- metadata$age
names(sex_vec) <- metadata$sample
names(age_vec) <- metadata$sample
tonsil_king$sex <- sex_vec[tonsil_king$Sample]
tonsil_king$age <- age_vec[tonsil_king$Sample]
tonsil_king$age_group <- "kid"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "Donor"] <- "donor_id"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "Sample"] <- "gem_id"
tonsil_king$hospital <- "Royal London"
tonsil_king$is_hashed <- "not_hashed"
tonsil_king$library_name <- tonsil_king$gem_id
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "Assay"] <- "assay"
tonsil_cnag$assay <- "3P"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "Status"] <- "status"
tonsil_cnag$status <- NA


# Calculate cell cycle scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
tonsil_cnag <- CellCycleScoring(
  tonsil_cnag,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)
tonsil_cnag$CC.Difference <- tonsil_cnag$S.Score - tonsil_cnag$G2M.Score


# Set IG columns to NA in our dataset
ig_columns <- str_subset(colnames(tonsil_king@meta.data), "^IG")
for (x in ig_columns) {
  tonsil_cnag@meta.data[[x]] <- NA
}
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "ChainStatus"] <- "chain_status"
tonsil_cnag$chain_status <- NA
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "ISOTYPE"] <- "isotype"
tonsil_cnag$isotype <- NA
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "SEQUENCE_ID"] <- "sequence_id"
tonsil_cnag$sequence_id <- NA
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "CLONE_ID"] <- "clone_id"
tonsil_cnag$clone_id <- NA


# Set annotation columns to NA in our dataset
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "CellType"] <- "cell_type"
tonsil_cnag$cell_type <- "unannotated"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "Lineage"] <- "lineage"
tonsil_cnag$lineage <- "unannotated"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "Subset"] <- "subset"
tonsil_cnag$subset <- "unannotated"
colnames(tonsil_king@meta.data)[colnames(tonsil_king@meta.data) == "MBC_Subset"] <- "MBC_subset"
tonsil_cnag$MBC_subset <- "unannotated"


# Merge
tonsil_king@meta.data <- tonsil_king@meta.data[, colnames(tonsil_cnag@meta.data)]
Idents(tonsil_cnag) <- "gem_id"
Idents(tonsil_king) <- "gem_id"
tonsil <- merge(x = tonsil_cnag, y = tonsil_king)


# Save
saveRDS(tonsil, path_to_save)
