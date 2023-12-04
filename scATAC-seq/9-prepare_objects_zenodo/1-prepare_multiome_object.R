# This script prepares the 10X Multiome object to share in Zenodo

# For ATAC, we need to:
#     1. Select all cells that we will include. That's a mix of the final object for the discovery cohort + newly sequenced multiome ATAC. We only include cells for which we have both ATAC + RNA info.
#     2. Define an annotation to use (e.g. UCSC)
#     3. Identify fragment files, only keep fragments from cells in the object. Save paths to fragment files to include in object.
#     4. Call peaks per cell type. The purpose is to identify peaks that are specific to all cell types, even rare ones. We might need to merge cell types to ensure we dont have too few cells
#     5. Create accessibility matrix based on the previous peaks (4) and the fragments file (3)
#     6. Normalize and perform dimensionality reduction


# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(glue)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(harmony)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) # for 100 Gb RAM


# Load object that contains RNA+multiome, subset to keep multiome only
seurat_rna <- readRDS(here("scRNA-seq/results/R_objects/final_clusters/20230911/20230911_tonsil_atlas_rna_seurat_obj.rds"))
seurat_rna <- seurat_rna[, seurat_rna$assay == "multiome"]


# Select cells
seurat_atac_all <- readRDS(here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_tonsil_atlas_atac_seurat_obj.rds"))
selected_cells <- c(
    seurat_rna$barcode[seurat_rna$barcode %in% seurat_atac_all$barcode], # these cells are from Multiome from the discovery cohort
    seurat_rna$barcode[seurat_rna$cohort_type == "validation"] # these cells were sequenced in the revision (Multiome, validation cohort)
)
seurat_rna <- seurat_rna[, selected_cells]


# Set fragments files
gem_ids <- unique(seurat_rna$gem_id)
fragment_objs <- map(gem_ids, \(gem_id) {
    path_frags <- glue("./fragments_multiome/{gem_id}_atac_fragments_with_prefix.tsv.gz")
    cells <- seurat_rna$barcode[seurat_rna$gem_id == gem_id]
    names(cells) <- cells
    frags_obj <- CreateFragmentObject(path = path_frags, cells = cells)
    frags_obj
})
names(fragment_objs) <- gem_ids


# Call peaks
peaks_list <- map(fragment_objs, \(frags_obj) {
  peaks <- CallPeaks(
    object = frags_obj,
    macs2.path = "~/anaconda3/envs/richter/bin/macs2"
  )
})
combined_peaks <- GenomicRanges::reduce(c(
  peaks_list[[1]],
  peaks_list[[2]],
  peaks_list[[3]],
  peaks_list[[4]],
  peaks_list[[5]],
  peaks_list[[6]],
  peaks_list[[7]],
  peaks_list[[8]],
  peaks_list[[9]],
  peaks_list[[10]],
  peaks_list[[11]],
  peaks_list[[12]]
))


# Create accessibility matrix
acc_mat <- FeatureMatrix(
  features = combined_peaks,
  fragment_objs,
  cells = seurat_rna$barcode
)
acc_mat <- acc_mat[, seurat_rna$barcode]


# Create Chromatin Assay and add to Seurat object 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
chrom_assay <- CreateChromatinAssay(
  counts = acc_mat,
  sep = c(":", "-"),
  fragments = fragment_objs,
  annotation = annotation
)
seurat_rna[["ATAC"]] <- chrom_assay


# Recall peaks by cell type to detect peaks enriched in rare populations,
# keep only standard chromosomes
DefaultAssay(seurat_rna) <- "ATAC"
peaks_cell_type <- CallPeaks(seurat_rna, group.by = "annotation_level_1")
peaks_cell_type <- keepStandardChromosomes(peaks_cell_type, pruning.mode = "coarse")
peaks_cell_type <- subsetByOverlaps(x = peaks_cell_type, ranges = blacklist_hg38_unified, invert = TRUE)


# Recreate accessibility matrix
acc_mat2 <- FeatureMatrix(fragments = fragment_objs, features = peaks_cell_type)
acc_mat2 <- acc_mat2[, colnames(seurat_rna)]
seurat_rna[["ATAC"]] <- CreateChromatinAssay(
  counts = acc_mat2,
  sep = c(":", "-"),
  fragments = fragment_objs,
  annotation = annotation
)

# Process ATAC-seq
DefaultAssay(seurat_rna) <- "ATAC"
seurat_rna <- seurat_rna %>%
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = "q0") %>%
  RunSVD()
seurat_rna <- RunHarmony(
  seurat_rna,
  reduction = "lsi",
  assay.use = "ATAC",
  reduction.save = "harmony_atac",
  project.dim = FALSE,
  group.by.vars = "donor_id"
)


# Process RNA-seq
DefaultAssay(seurat_rna) <- "RNA"
seurat_rna <- seurat_rna %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(
    reduction = "pca",
    group.by.vars = "donor_id",
    reduction.save = "harmony_rna",
    assay.use = "RNA",
    project.dim = FALSE
  )


# Multimodal nearest neighbors and UMAP
seurat_rna <- FindMultiModalNeighbors(
  seurat_rna,
  reduction.list = list("harmony_rna", "harmony_atac"),
  dims.list = list(1:30, 2:30)
)
seurat_rna <- RunUMAP(
  seurat_rna,
  nn.name = "weighted.nn",
  verbose = TRUE
)


# Calculate ATAC QC variables
DefaultAssay(seurat_rna) <- "ATAC"
seurat_rna <- NucleosomeSignal(object = seurat_rna)
seurat_rna <- TSSEnrichment(object = seurat_rna, fast = FALSE)
DefaultAssay(seurat_rna) <- "RNA"


# Remove clutter
seurat_rna@reductions$harmony <- NULL
seurat_rna$age_group[seurat_rna$age_group == "kid"] <- "child"


# Correct annotation
seurat_rna$annotation_20230508[seurat_rna$annotation_level_1 == "FDC"] <- "FDC"
seurat_rna$annotation_20230508[seurat_rna$annotation_level_1 == "epithelial"] <- "epithelial"
seurat_rna$annotation_20230508[seurat_rna$annotation_level_1 == "PDC"] <- "PDC"


# Correct metadata
seurat_rna$cause_for_tonsillectomy[seurat_rna$cause_for_tonsillectomy == "spleep_apnea"] <- "sleep apnea"
seurat_rna$cause_for_tonsillectomy[seurat_rna$cause_for_tonsillectomy == "sleep_apnea"] <- "sleep apnea"

# Save
saveRDS(seurat_rna, here("scRNA-seq/results/R_objects/final_clusters/20230911/20230911_tonsil_atlas_multiome_seurat_obj.rds"))