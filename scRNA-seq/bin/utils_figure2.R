# This script contains parameters, variables and functions used to create
# figure 2


# PATHS
path_save_umap_fig2 <- here("results/paper/figures/figure_2_main_umap.pdf")
path_save_heatmap_rna <- here("results/paper/figures/figure_2_heatmap_rna.pdf")
path_save_heatmap_rna_leg <- here("results/paper/figures/figure_2_heatmap_rna_legend.pdf")
path_save_umaps_act <- here("results/paper/figures/figure_2_umaps_tf_activities.pdf")
path_save_umaps_act_leg <- here("results/paper/figures/figure_2_umaps_tf_activities_legend.pdf")
path_to_save_figure_2_scTCR_seq <- here("results/paper/figures/figure_2_tcr.pdf")
path_save_heatmaps_tfh_rna <- here("results/paper/figures/figure_2_heatmaps_tfh_rna.pdf")
path_save_heatmaps_tfh_atac <- here("results/paper/figures/figure_2_heatmaps_tfh_atac.pdf")
path_save_heatmaps_tfh_atac_legend <- here("results/paper/figures/figure_2_heatmaps_tfh_atac_legend.pdf")
path_save_coverage_plot <- here("results/paper/figures/figure_2_coverage_plot_tfh.pdf")
path_save_umaps_acc <- here("results/paper/figures/figure_2_umaps_bcl_accessibilty.pdf")
path_save_umaps_acc_leg <- here("results/paper/figures/figure_2_umaps_bcl_accessibilty_legend.pdf")
path_to_save_th_main <- here("results/paper/figures/figure_2_th_nebulosa.pdf")
path_save_umaps_cite <- here("results/paper/figures/figure_2_umaps_cite.pdf")
path_save_umaps_cite_leg <- here("results/paper/figures/figure_2_umaps_cite_legend.pdf")


# FUNCTIONS
downsample_cells2 <- function(seurat_obj,
                              n_cells_total,
                              var = "annotation_paper") {
  cells <- colnames(seurat_obj)
  categories <- unique(seurat_obj@meta.data[[var]])
  n_cells <- round(n_cells_total / length(categories))
  selected_cells <- purrr::map(categories, function(x) {
    sampled <- sample(
      cells[seurat_obj$annotation_paper == x],
      size = n_cells,
      replace = FALSE
    )
    sampled
  })
  selected_cells <- unlist(selected_cells)
  names(selected_cells) <- NULL
  seurat_sub <- subset(seurat_obj, cells = selected_cells)
  seurat_sub
}

correlation_heatmap2 <- function(
  se,
  genes = NULL,
  feats = NULL,
  assay = "Spatial",
  slot = "data",
  cor_method = "pearson",
  text_size = 5,
  type = "upper",
  ...) {
  
  if (!is.null(genes)) {
    
    # Extract expression matrix
    expr_mtrx <- Seurat::GetAssayData(
      object = se,
      assay = assay,
      slot = slot)
    
    # Return warning for those genes not found in the expression matrix
    if (sum(genes %in% rownames(expr_mtrx)) != length(genes)) {
      gene_out <- genes[! genes %in% rownames(expr_mtrx)]
      gene_out <- paste(gene_out, collapse = ", ")
      warning(glue::glue("{gene_out} not found in assay {assay}, slot {slot}"))
    }
    
    # Subset expression matrix to genes of interest
    genes_sub <- genes[genes %in% rownames(expr_mtrx)]
    
    # Deal with behaviour of creating a matrix with just 1 row
    if (length(genes_sub) == 1) {
      mtrx_genes <- as.matrix(expr_mtrx[genes_sub, ])
      colnames(mtrx_genes) <- genes_sub
    } else if (length(genes_sub) > 1) {
      mtrx_genes <- t(as.matrix(expr_mtrx[genes_sub, ]))
    }
  }
  
  # Extract features of interest
  if (!is.null(feats)) {
    mtrx_feats <- se@meta.data[, feats]
  }
  
  # Sanity check and combining genes and feats into mtrx if both are defined
  if (is.null(genes) & is.null(feats)) {
    
    stop("Need to pass either genes or metadata features")
    
  } else if ((! is.null(genes)) & is.null(feats)) {
    
    mtrx <- mtrx_genes
    
  } else if (is.null(genes) & (! is.null(feats))) {
    
    mtrx <- mtrx_feats
    
  } else {
    
    mtrx <- cbind(mtrx_genes, mtrx_feats)
    
  }
  
  # Remove features that are all 0
  mtrx <- mtrx[, colSums(mtrx) > 0]
  mtrx_cor <- cor(as.matrix(mtrx))
  
  # Compute correlation P-value
  p_mat <- corrplot::cor.mtest(mat = mtrx,
                               conf_int = 0.95,
                               method = cor_method)
  
  # Add new line in names over 30 characters long
  colnames(mtrx_cor) <- stringr::str_wrap(string = colnames(mtrx_cor),
                                          width = 30)
  rownames(mtrx_cor) <- stringr::str_wrap(string = rownames(mtrx_cor),
                                          width = 30)
  
  # Plot correlation matrix as a heatmap
  ggcorrplot::ggcorrplot(
    corr = mtrx_cor,
    p.mat = p_mat[[1]],
    hc.order = TRUE,
    type = type,
    insig = "blank",
    lab = FALSE,
    outline.col = "lightgrey",
    method = "square",
    colors = c("#6D9EC1", "white", "#E46726"),
    legend.title = glue::glue("Correlation\n({cor_method})"),
    ...) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = text_size),
      axis.text.y = ggplot2::element_text(vjust = 0.5, size = text_size))
}
theme_nothing2 <- function() {
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
}
no_legend <- function() {
  theme(legend.position = "none")
}


# PATHS
path_to_auc_mtx <- here("data/raw_data_figures/all_Tcells_auc.csv")
path_to_mat_rna_bcl6 <- here("data/raw_data_figures/scATAC_CD4_T/files_plots/matrix_RNA_genes.tsv")
path_to_mat_atac_bcl6 <- here("data/raw_data_figures/scATAC_CD4_T/files_plots/matrix_ATAC_genes.tsv")


# path_to_cite <- "~/Google Drive/My Drive/single_cell/PhD/B_cell_atlas/tonsil_atlas/tonsil_atlas_annotation/current/preprint/CD4_T/CD4_tonsil_cite_seq_annotated_tcr_analysis_updated_obj.rds"
# path_to_mat_rna <- "/home/rmassonix/Desktop/figure2_ATAC/matrix_RNA_genes.tsv"
# path_to_mat_atac <- "/home/rmassonix/Desktop/figure2_ATAC/matrix_ATAC_genes.tsv"
# path_to_atac <- "~/Desktop/figure2_ATAC/CD4T_integration_peak_calling_level_5.rds"
# path_to_th <- "~/Google Drive/My Drive/single_cell/PhD/B_cell_atlas/tonsil_atlas/tonsil_atlas_annotation/current/preprint/CD4_T/CD4_T_Th_level_5_annotated.rds"
# path_to_auc_mtx <-  "~/Desktop/tonsil_atlas/scRNA-seq/results/all_Tcells_auc.csv"
# path_to_spata <- "~/Google Drive/My Drive/single_cell/PhD/B_cell_atlas/tonsil_atlas/tonsil_atlas_annotation/current/ST/spata-esvq52_nluss5.rds"
# path_to_spatial <- "~/Google Drive/My Drive/single_cell/PhD/B_cell_atlas/tonsil_atlas/tonsil_atlas_annotation/current/ST/se_sub.rds"
# path_to_save <- "~/Desktop/tonsil_atlas/scRNA-seq/results/plots/tonsil_atlas_figure2.pdf"


# COLORS
colors_rna <- c(
  "Naive" = "gray88",
  "CM Pre-non-Tfh" = "gray71",
  "CM PreTfh" = "gray59",
  "Tfh T:B border" = "#67a9cf",
  "Tfh-LZ-GC" = "#3690c0",
  "GC-Tfh-SAP" = "#02818a",
  "GC-Tfh-0X40"= "#016c59",
  "Tfh-Mem" = "#014636",
  "T-Trans-Mem" = "#fd8d3c",
  "T-Eff-Mem" = "#e31a1c",
  "T-helper" = "#800026",
  "Eff-Tregs" = "#df65b0",
  "Eff-Tregs-IL32" = "#e7298a",
  "Tfr" = "#ce1256"
)
colors_th <- c(
  "Th0" = "#7fc97f",
  "Th1" = "#beaed4",
  "Th2" = "#fdc086",
  "Th17" = "#ffff99",
  "Th1/Th17" = "#386cb0",
  "Th22" = "#f0027f"
)
hm_colors <- viridis::inferno(n = 100)
colors_tcr <- c("#c3c3c3", "#f0b635")


# GENES OF INTEREST (GOI)
goi_rna <- list(
  "Naive" = c("BACH2", "LEF1", "CCR7", "NOSIP"),
  "CM Pre-non-Tfh" = c("ANXA1", "IL7R", "ITGB1", "S100A10", "S100A4"),
  "CM PreTfh" = c("ANK3", "ZBTB16", "TXK", "PCNX1", "FKBP5"),
  "Tfh T:B border" = c("TOX", "PTPN13", "SLC9A9", "IL6ST"),
  "Tfh-LZ-GC" = c("PASK", "ST8SIA1", "IFITM1", "TOX2", "ID3"),
  "GC-Tfh-SAP" = c("SH2D1A", "CXCL13", "PDCD1", "POU2AF1"),
  "GC-Tfh-0X40"= c("CD200", "EGR2", "IL21", "TNFRSF4", "BATF", "IRF4"),
  "Tfh-Mem" = c("CEP128", "MAF", "IKZF3", "TNFRSF1B"),
  "T-Trans-Mem" = c("KLRB1", "RORA", "RUNX2", "CCR6"),
  "T-Eff-Mem" = c("HECW2", "SESN3", "PTMS", "EGLN3"),
  "T-helper" = c("LAG3", "IL17A", "CCL4", "IL26", "CCL20", "CXCR6"),
  "Eff-Tregs" = c("IL2RA", "CYTOR", "RGS1"),
  "non-GC-Tf-regs" = c("IKZF2", "FOXP3", "IL32"),
  "GC-Tf-regs" = c("FCRL3", "CLNK", "SESN3")
)
goi_rna_supp <- c("ICOS", "CXCR5", "TCF7", "IL12RB2", "TGFBR1", "CD40LG",
                  "BTLA", "TNFRSF4", "IL4")
# goi_rna <- list(
#   "Naive" = c("CCR7", "LEF1", "NOSIP"),
#   "CM Pre-non-Tfh" = c("ANXA1", "IL7R", "S100A4"),
#   "CM PreTfh" = c("TXK"),
#   "Tfh T:B border" = c("ICOS", "PDCD1", "CXCR5", "CXCL13", "TOX", "IL6ST"),
#   "Tfh-LZ-GC" = c("IFITM1", "TOX2"),
#   "GC-Tfh-SAP" = c("SH2D1A"),
#   "GC-Tfh-0X40"= c("CD200", "EGR2", "TNFRSF4"),
#   "Tfh-Mem" = c("CEP128", "KLRB1"),
#   "T-Trans-Mem" = c("RORA", "CCR6"),
#   "T-Eff-Mem" = c("SESN3"),
#   "T-helper" = c("LAG3", "IL17A"),
#   "Tregs" = c("CTLA4", "FOXP3", "IL2RA", "IKZF2")
#   # "non-GC-Tf-regs" = c("IKZF2", "FOXP3", "IL32"),
#   # "GC-Tf-regs" = c("FCRL3", "CLNK", "SESN3")
# )


goi_cite <- c("CD45RA", "CD62L", "CD127-(IL-7Ralpha)", "CD45RO",
              "CD279-(PD-1)", "CD278-(ICOS)", "CD185-(CXCR5)", "CD161",
              "CD25")
# goi_cite <- c("CD28", "CD27", "CD197-(CCR7)", "CD62L", "CD45RA", "CD45RO",
#               "CD95-(Fas)", "CD127-(IL-7Ralpha)", "CD122-(IL-2Rbeta)", "CD38",
#               "CD58-(LFA-3)", "CD185-(CXCR5)", "CD194-(CCR4)", "CD196-(CCR6)",
#               "CD183-(CXCR3)", "CD161", "CD25", "CD278-(ICOS)", "CD252-(OX40L)")
# goi_cite <- c(
#   naive_central_memory = c("CD49f", "CD62L", "CD197-(CCR7)", "CD127-(IL-7Ralpha)", "CD45RA"),
#   tfh = c("CD45RO", "CD279-(PD-1)", "CD278-(ICOS)", "TIGIT-(VSTM3)"),
#   teff = c("CD161", "CD28", "CD95-(Fas)", "CD39"),
#   treg = c("CD25", "CD71")
# )


goi_th_main <- c("IL7R", "IFNG", "CXCR3", "CCR4", "IL17F", "IL17A",
                 "IL22", "IL10")
goi_th_supp <- c("GATA3", "TBX21", "GZMK", "CCR5", "IL17RB", "IL4", "CCR6",
                 "RORC", "IL26", "TNF")

goi_spatial <- list(
  "DC" <- c("CCR7", "HLA-DRA", "CLEC9A"),
  "B-cognate" <- c("CD40", "CD80", "CD86", "BCL6", "ICOSLG"),
  "CD4-naive" <- c("CCR7", "CD4", "SELPLG", "TCF7", "LEF1", "STAT1", "STAT3"),
  "CD4-Tfh" <- c("BCL6", "CD40LG", "CXCR5", "IL6R", "ICOS", "IL21", "CD28", "SELPLG", "GPR183"),
  "CD4-non-Tfh" <- c("PRDM1", "IL2RA"),
  "CD4-GC-Tfh" <- c("SH2D1A", "CD200", "BTLA", "PDCD1", "CXCR5", "GPR183", "IL4"),
  "Chemokines" <- c("IL2", "IL6", "IL12A", "IL21", "IL10"),
  "FDC" <- c("CR1", "CR2", "CXCL13"),
  "Proliferation" <- c("MKI67", "CDK1", "TOP2A"),
  "Lymph Endothelial" <- c("PECAM1", "CLDN5", "PROX1", "LYVE1"),
  "T cells" <- c("CD3D", "CD3E", "CD8A", "CD8B"),
  "B cells" <- c("MS4A1", "CD79A", "CD40", "ICAM1", "ICOSLG", "CD80"),
  "DZ" <- c("OAZ1", "AICDA", "H3", "MKI67", "POLH"),
  "LZ" <- c("LAG3", "ITGB8", "PDCD1", "TIGIT", "BCL2", "PECAM1", "LY6E",
            "CD276", "HLA-DRB1", "PSMB10", "TNF", "ARG1", "HLA-E", "STAT1"),
  "naDC" <- c("ITGAX", "XCR1")
)
