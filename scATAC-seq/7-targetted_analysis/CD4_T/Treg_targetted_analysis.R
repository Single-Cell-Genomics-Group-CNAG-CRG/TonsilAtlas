# Load packages
library(Seurat)
library(Signac)
library(scCustomize)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
library(patchwork)
library(tidyverse)
library(here)


# Souce utilities
source("scRNA-seq/bin/utils_final_clusters.R")
source("scRNA-seq/bin/utils_figure2.R")

# Read data
treg_levels <- c("Eff-Tregs", "non-GC-Tf-regs", "GC-Tf-regs")
seurat_rna <- readRDS(path_to_save_cd4)
seurat_rna <- subset(seurat_rna, idents = treg_levels)
seurat <- readRDS(path_to_save_atac_cd4)
seurat <- subset(seurat, idents = treg_levels)


# Call peaks
peaks_level3 <- CallPeaks(
  object = seurat,
  group.by = "annotation_20220215",
  macs2.path = "/home/rmassonix/miniconda3/envs/macs2/bin/macs2"
)
peaks_level3 <- keepStandardChromosomes(peaks_level3, pruning.mode = "coarse")
peaks_level3 <- subsetByOverlaps(
  peaks_level3,
  ranges = blacklist_hg38_unified,
  invert = TRUE
)


# Call accessibility matrix
annotation <- GetGRangesFromEnsDb(
  ensdb = EnsDb.Hsapiens.v86,
  standard.chromosomes = TRUE
)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seurat),
  features = peaks_level3,
  cells = colnames(seurat)
)
seurat[["peaks_redefined"]] <- CreateChromatinAssay(
  counts = macs2_counts, 
  genome = "hg38",
  fragments = Fragments(seurat),
  annotation = annotation
)
DefaultAssay(seurat) <- "peaks_redefined"


# Dimensionality reduction
seurat@assays$peaks_level_5 <- NULL
seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = "q50")
seurat <- RunSVD(seurat)
DepthCor(seurat)
seurat <- RunHarmony(
  seurat,
  group.by.vars = "assay",
  dims = 2:20,
  reduction = "lsi",
  assay.use = "peaks_redefined"
)
seurat <- RunUMAP(object = seurat, reduction = "harmony", dims = 2:20)
DimPlot(seurat, group.by = "annotation_20220215")


# Motif analysis
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species = "Homo sapiens", all_versions = FALSE)
)
seurat <- AddMotifs(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
seurat <- RunChromVAR(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
DefaultAssay(seurat) <- "chromvar"
saveRDS(
  seurat,
  here("scATAC-seq/results/R_objects/Targetted_analysis/CD4_T/20220412_seurat_treg_atac.rds")
)
# seurat <- readRDS(here("scATAC-seq/results/R_objects/Targetted_analysis/CD4_T/20220412_seurat_treg_atac.rds"))
differential_activity <- FindAllMarkers(
  object = seurat,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
differential_activity$motif_name <- unlist(seurat@assays$peaks_redefined@motifs@motif.names[rownames(differential_activity)])
differential_activity$motif <- rownames(differential_activity)

# a <- VlnPlot(seurat, features = "MA1508.1", group.by = "annotation_20220215", pt.size = 0)
# a <- a +
#   labs(title = "IKZF1", y = "Accessibility") +
#   theme(legend.position = "none", axis.title.x = element_blank())
# b <- MotifPlot(
#   object = seurat,
#   motifs = "MA1508.1",
#   assay = 'peaks_redefined'
# ) + theme(plot.title = element_blank())
# a | b


# Markers RNA
differential_expression <- FindAllMarkers(
  seurat_rna,
  logfc.threshold = 0.75,
  only.pos = TRUE
)
overlapping <- intersect(differential_activity$motif_name, differential_expression$gene)
overlapping <- overlapping[overlapping != "BACH1"] # BACH1 is in different cell types
vln_plots_expr <- Stacked_VlnPlot(
  seurat_rna,
  features = overlapping,
  group.by = "annotation_20220215",
  colors_use = colors_rna[treg_levels],
  plot_spacing = 0.05
)
overlapping_motifs <- differential_activity[differential_activity$motif_name %in% overlapping, "motif"]
vln_plots_acc <- Stacked_VlnPlot(
  seurat,
  features = overlapping_motifs,
  group.by = "annotation_20220215",
  colors_use = colors_rna[treg_levels],
  plot_spacing = 0.05
)
motif_plots <- map(overlapping_motifs, \(. ) MotifPlot(seurat, ., assay = "peaks_redefined"))
motif_plots_arr <- wrap_plots(motif_plots, ncol = 1)
motif_plots_arr | vln_plots_acc | vln_plots_expr
MotifPlot(object = seurat, overlapping_motifs, assay = "peaks_redefined")



# # Transcription factor footprinting
# seurat <- Footprint(
#   object = seurat,
#   motif.name = c("MA0495.3"),
#   genome = BSgenome.Hsapiens.UCSC.hg38
# )
# 
# # plot the footprint data for each group of cells