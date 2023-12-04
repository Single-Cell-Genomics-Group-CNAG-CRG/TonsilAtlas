# This script plots all the results for the figure 7

# Load packages
library(Seurat)
library(tidyverse)
library(ggrastr)
library(ggpubr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(here)
library(glue)


# Source functions and colors
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure2.R"))
source(here("scRNA-seq/bin/utils_figure3.R"))
goi_cytotoxic <- goi_rna
source(here("scRNA-seq/bin/utils_figure4.R"))

# Read data
cell_types <- c("NBC_MBC", "GCBC", "pc", "Cytotoxic", "myeloid")
seurat_list <- map(cell_types, \(x) {
  file <- list.files(
    here("scRNA-seq/results/R_objects/seurat_objects_revision"),
    pattern = glue("^5.*{x}"),
    full.names = TRUE
  )
  print(file)
  seurat_obj <- readRDS(file)
  seurat_obj
})
names(seurat_list) <- cell_types
names(seurat_list)[names(seurat_list) == "pc"] <- "PC"
cell_types[cell_types == "pc"] <- "PC"


# Plot UMAPs
umaps <- map(cell_types, \(x) {
  seurat_obj <- seurat_list[[x]]
  colors <- colors_20230508[[x]]
  seurat_obj$annotation_20230508 <- droplevels(seurat_obj$annotation_20230508)
  p <- Embeddings(seurat_obj, "umap") %>%
    as.data.frame() %>%
    mutate(
      annotation = seurat_obj$annotation_20230508,
      cohort_type = seurat_obj$cohort_type
    ) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = annotation))
  if (x %in% c("NBC_MBC", "GCBC", "Cytotoxic", "PC")) {
    p <- p + geom_point(shape = ".", alpha = 0.9)
  } else if (x == "myeloid") {
    p <- p + geom_point(size = 0.05)
  }
  p <- p +
    facet_wrap(~cohort_type) +
    scale_color_manual(values = colors, breaks = names(colors)) +
    theme_classic() +
    theme(
      strip.text = element_text(size = 6),
      strip.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.position = "bottom",
      legend.spacing = unit(0, "cm"),
      legend.box.spacing = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.key.height = unit(0, "cm"),
      legend.key.width = unit(0, "cm"),
      plot.title = element_text(hjust = 0.5, size = 7),
      axis.title = element_text(size = 6),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    ) +
    guides(color = guide_legend(nrow = 4, override.aes = list(shape = 16, size = 2)))
  p <- rasterize(p, dpi = 600)
  p
})
names(umaps) <- cell_types
for (cell_type in cell_types) {
  print(cell_type)
  ggsave(
    filename = here(glue("results/paper/figures/supplementary_figure_umap_{cell_type}_annotation_validation_cohort.pdf")),
    plot = umaps[[cell_type]],
    width = 10,
    height = 6,
    units = "cm"
  )
}

for (cell_type in cell_types) {
  print(cell_type)
  legend <- as_ggplot(get_legend(umaps[[cell_type]]))
  ggsave(
    filename = here(glue("results/paper/figures/supplementary_figure_umap_{cell_type}_annotation_validation_cohort_legend.pdf")),
    plot = legend,
    width = 18,
    height = 18,
    units = "cm"
  )
}


# Heatmap markers
my_palette <- colorRamp2(
  c(0, 0.25, 0.5, 0.75, 1),
  c("#4976b3", "#a7d5e6", "#f3ffc7", "#ffc685", "#f20034")
)
in2mm <- 25.4
goi_nbc_mbc <- c("TCL1A", "FCER2", "CD27", "TNFRSF13B", "IGHD", "IGHM", "IGHA1",
                "IGHG1", "CD69", "EGR2", "CCL3", "MYC", "BATF", "CCND2", "NFKB1",
                "NFKB2", "RELB", "CD40", "BCL2A1", "MIR155HG", "EBI3", "TRAF1",
                "MEF2B", "RGS13", "MKI67", "TOP2A", "LY9", "LILRA4", "IFIT1", 
                "IFIT3", "FCRL4", "FCRL5")
goi_gcbc <- c("CXCR4", "AICDA", "FOXP1", "MME", "CD83", "LMO2", "POLA1", "HIST1H4C",
              "MKI67", "TOP2A", "CDC20", "CCNB1", "STMN1", "HMGB2", "CD40", "BCL2A1",
              "MIR155HG", "EBI3", "TRAF1", "NFKB1", "NFKB2", "RELB", "MYC", "BATF",
              "CCR6", "CELF2", "BANK1", "PRDM1", "XBP1")
goi_pc <- c("SUGCT", "AICDA", "CXCR4", "LMO2", "CD83", "BCL2A1", "BCL6", "IRF8",
            "MEF2B", "MS4A1", "PAX5", "PRDM1", "XBP1", "IRF4", "SLAMF7", "SSR4",
            "MZB1", "DERL3", "CREB3L2", "FKBP11", "IGHG1", "IGHG2", "IGHG3", "IGHG4",
            "IGHA1", "IGHA2", "IGHM", "IGHD", "CCDC50", "FCMR", "CD9", "CD44", "H2AFZ",
            "MKI67", "TUBA1B", "TOP2A")
goi_myeloid <- c(goi_myeloid, "SELENOP", "APOE", "MMP9", "MMP12", "C1QA", "ITGAX", "ZEB2")
goi <- list(
  "NBC_MBC" = goi_nbc_mbc,
  "GCBC" = goi_gcbc,
  "PC" = goi_pc,
  "Cytotoxic" = goi_cytotoxic,
  "myeloid" = goi_myeloid

)
for (cell_type in cell_types) {
  print(cell_type)
  seurat_obj <- seurat_list[[cell_type]][, seurat_list[[cell_type]]$cohort_type == "validation"]
  seurat_obj$annotation_20230508 <- droplevels(seurat_obj$annotation_20230508)
  colors <- colors_20230508[[cell_type]]
  avgexpr_mat <- AverageExpression(
    seurat_obj,
    features = goi[[cell_type]],
    assays = "RNA",
    return.seurat = FALSE,
    group.by = "annotation_20230508",
    slot = "data"
  )$RNA
  input_mat <- apply(avgexpr_mat, 1, function(x) (x - min(x)) / diff(range(x)))
  probs_l <- map(levels(seurat_obj$annotation_20230508), \(x) {
    x <- seurat_obj$annotation_20230508_probability[seurat_obj$annotation_20230508 == x & seurat_obj$cohort_type == "validation"]
    x
  })
  names(probs_l) <- levels(seurat_obj$annotation_20230508)
  anno_row <- rowAnnotation(
    cell_type = levels(seurat_obj$annotation_20230508),
    annotation_probability = anno_boxplot(
      x = probs_l,
      width = unit(1, "cm"),
      size = unit(1, "pt")
    ),
    annotation_name_gp = gpar(fontsize = 6) ,
    col = list("cell_type" = colors),
    show_legend = c(FALSE, FALSE, FALSE)
  )
  path_save_heatmap_rna <- here(glue("results/paper/figures/supplementary_figure_heatmap_{cell_type}_markers_validation_cohort.pdf"))
  pdf(path_save_heatmap_rna, width = (110 / in2mm), height = (55 / in2mm), paper = "special")
  heatmap <- Heatmap(
    input_mat,
    cluster_rows = FALSE,
    show_heatmap_legend = FALSE,
    cluster_columns = FALSE,
    column_names_gp = gpar(fontsize = 5.5),
    col = my_palette,
    name = "Expression",
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 7),
    right_annotation = anno_row,
    heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 6))
  )
  draw(heatmap)
  dev.off()
}



