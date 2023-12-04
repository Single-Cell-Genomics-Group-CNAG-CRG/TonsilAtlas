# This script plots all the results for the figure 7

# Load packages
library(Seurat)
library(tidyverse)
library(ggrastr)
library(here)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# Source functions and colors
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure2.R"))


# Read data
seurat <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/7-20230529_seurat_rna_discovery_validation_cohorts.rds"))
cd4 <- readRDS(here("scRNA-seq/results/R_objects/seurat_objects_revision/5.8-CD4_seurat_object_discovery_validation_cohorts_annotation.rds"))


# Plot UMAP level 1 split by cohort type
seurat$annotation_level_1 <- factor(
  seurat$annotation_level_1,
  names(colors_20230508$level_1)
)
umap_level_1 <- Embeddings(seurat, "umap") %>%
  as.data.frame() %>%
  mutate(
    annotation_level_1 = seurat$annotation_level_1,
    cohort_type = seurat$cohort_type
  ) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = annotation_level_1)) +
  geom_point(shape = ".", alpha = 0.9) +
  facet_wrap(~cohort_type) +
  scale_color_manual(values = colors_20230508$level_1, breaks = names(colors_20230508$level_1)) +
  # labs(title = str_c(umap_title, "cells", sep = " "), x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  theme(
    strip.text = element_text(size = 6),
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
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2)))

umap_level_1 <- rasterize(umap_level_1, dpi = 600)
ggsave(
  filename = here("results/paper/figures/figure_7_umap_annotation_level_1_split_by_cohort.pdf"),
  plot = umap_level_1,
  width = 14,
  height = 7,
  units = "cm"
)
Idents(seurat) <- "annotation_level_1"
# markers_level_1 <- FindAllMarkers(
#   seurat[, seurat$cohort_type == "discovery"],
#   only.pos = TRUE,
#   logfc.threshold = 0.75
# )
markers_level_1 <- c(
  "BANK1", "IGHD",
  "RGS13", "MEF2B", "MYBL1",
  "JCHAIN", "MZB1", "XBP1",
  "TRAC", "CD3E", "IL7R",
  "GZMK", "GZMA", "NKG7",
  "LYZ", "CST3", "APOE",
  "FDCSP", "CLU", "CXCL13",
  "LILRA4", "IL3RA", "ITM2C",
  "SPRR3", "KRT19", "KRT13",
  "VPREB1", "RAG1"
)
seurat <- seurat[, seurat$cohort_type == "validation"]
avgexpr_mat <- AverageExpression(
  seurat[, seurat$cohort_type == "validation"],
  features = markers_level_1,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "annotation_level_1",
  slot = "data"
)$RNA
input_mat <- apply(avgexpr_mat, 1, function(x) (x - min(x)) / diff(range(x)))
cell_types <- levels(seurat$annotation_level_1)
cell_types <- cell_types[cell_types != "preTC"]
probs_l <- map(cell_types, \(x) {
  x <- seurat$annotation_level_1_probability[seurat$annotation_level_1 == x & seurat$cohort_type == "validation"]
  x
})
names(probs_l) <- cell_types
total_cells_cell_type <- table(seurat$annotation_level_1)
total_cells_cell_type <- total_cells_cell_type[names(total_cells_cell_type) != "preTC"]
anno_row <- rowAnnotation(
  cell_type = cell_types,
  annotation_probability = anno_boxplot(
    x = probs_l,
    width = unit(1, "cm"),
    size = unit(1, "pt")
  ),
  number_of_cells = anno_barplot(
    as.matrix(total_cells_cell_type),
    add_numbers = FALSE,
    ylim = c(0, 100000),
    border = FALSE,
    width = unit(1, "cm"),
    numbers_gp = gpar(fontsize = 6),
    annotation_name_gp = gpar(fontsize = 6),
    show_legend = FALSE
  ),
  annotation_name_gp = gpar(fontsize = 6) ,
  col = list("cell_type" = colors_20230508$level_1)
)
my_palette <- colorRamp2(
  c(0, 0.25, 0.5, 0.75, 1),
  c("#4976b3", "#a7d5e6", "#f3ffc7", "#ffc685", "#f20034")
)
in2mm <- 25.4
path_save_heatmap_rna <- here("results/paper/figures/figure_7_heatmap_level_1_validation_cohort.pdf")
pdf(path_save_heatmap_rna, width = (140 / in2mm), height = (55 / in2mm), paper = "special")
Heatmap(
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
dev.off()




# CD4 T cells
umap_cd4 <- Embeddings(cd4, "umap") %>%
  as.data.frame() %>%
  mutate(
    annotation_20230508 = droplevels(cd4$annotation_20230508),
    cohort_type = cd4$cohort_type
  ) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = annotation_20230508)) +
  geom_point(shape = ".", alpha = 0.9) +
  facet_wrap(~cohort_type) +
  scale_color_manual(values = colors_20230508$CD4, breaks = names(colors_20230508$CD4)) +
  # labs(title = str_c(umap_title, "cells", sep = " "), x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  theme(
    strip.text = element_text(size = 6),
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
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2)))
umap_cd4 <- rasterize(umap_cd4, dpi = 600)
ggsave(
  filename = here("results/paper/figures/figure_7_umap_annotation_CD4_T_split_by_cohort.pdf"),
  plot = umap_cd4,
  width = 14,
  height = 7,
  units = "cm"
)


cd4 <- cd4[, cd4$cohort_type == "validation"]
goi_rna <- unlist(goi_rna)
names(goi_rna) <- NULL
goi_rna <- c(goi_rna, "TOP2A", "MKI67", "PCNA")
goi_rna <- goi_rna[!(goi_rna %in% c("LEF1", "S100A10", "TXK", "ID3", "TNFRSF4", "IL26", "CCL4"))]
avgexpr_mat <- AverageExpression(
  cd4,
  features = goi_rna,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "annotation_20230508",
  slot = "data"
)$RNA
input_mat <- apply(avgexpr_mat, 1, function(x) (x - min(x)) / diff(range(x)))
cell_types <- levels(cd4$annotation_20230508)
cell_types <- cell_types[cell_types != "unannotated"]
probs_l <- map(cell_types, \(x) {
  x <- cd4$annotation_20230508_probability[cd4$annotation_20230508 == x & cd4$cohort_type == "validation"]
  x
})
names(probs_l) <- cell_types
total_cells_cell_type <- table(cd4$annotation_20230508)
total_cells_cell_type <- total_cells_cell_type[names(total_cells_cell_type) != "unannotated"]
anno_row <- rowAnnotation(
  cell_type = cell_types,
  annotation_probability = anno_boxplot(
    x = probs_l,
    width = unit(1, "cm"),
    size = unit(1, "pt")
  ),
  number_of_cells = anno_barplot(
    as.matrix(total_cells_cell_type),
    add_numbers = FALSE,
    ylim = c(0, 18000),
    border = FALSE,
    width = unit(1, "cm"),
    numbers_gp = gpar(fontsize = 6),
    annotation_name_gp = gpar(fontsize = 6),
    show_legend = FALSE
  ),
  annotation_name_gp = gpar(fontsize = 6) ,
  col = list("cell_type" = colors_20230508$CD4)
)
my_palette <- colorRamp2(
  c(0, 0.25, 0.5, 0.75, 1),
  c("#4976b3", "#a7d5e6", "#f3ffc7", "#ffc685", "#f20034")
)
in2mm <- 25.4
path_save_heatmap_rna <- here("results/paper/figures/figure_7_heatmap_CD4_T_validation_cohort.pdf")
pdf(path_save_heatmap_rna, width = (150 / in2mm), height = (62.5 / in2mm), paper = "special")
Heatmap(
  input_mat,
  cluster_rows = FALSE,
  show_heatmap_legend = FALSE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 5),
  col = my_palette,
  name = "Expression",
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 7),
  right_annotation = anno_row,
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 6))
)
dev.off()



