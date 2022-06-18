# This script plots all panels related with myeloid cells for figure 4

# Load packages
library(Seurat)
library(UCell)
library(here)
library(ggpubr)
library(tidyverse)


# Source utils
source(here("scRNA-seq/bin/utils_final_clusters.R"))
source(here("scRNA-seq/bin/utils_figure4.R"))


# Read data
seurat <- readRDS(path_to_save_myeloid)
seurat$annotation_paper <- factor(
  seurat$annotation_20220215,
  levels = names(colors_myel)
)



# Extended dot plot DC
Idents(seurat) <- "annotation_paper"
goi_dc <- c("HLA-DPA1", "CADM1", "CAMK2D", "IDO1", "CLNK",
            "FCER1A", "CLEC10A", "HLA-DQB1", "HLA-DPB1", "S100A8",
            "S100A9", "ANXA1", "FTL", "SERPINA1", "LILRA4", "LST1", "TCF7L2",
            "AXL", "PPP1R14A", "CD22", "DAB2", "IL1RN",
            "IL4I1", "CD83", "AREG", "IFI30", "IL1B")
dc_levels <- c("DC1 precursor", "DC1 mature", "DC2", "DC3",
               "DC4", "DC5", "IL7R DC")
dc <- seurat[, seurat$annotation_20220215 %in% dc_levels]
dot_plot_dc <- DotPlot(dc, features = rev(goi_dc)) +
  coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1))
legend_dot_plot <- as_ggplot(get_legend(dot_plot_dc))
dot_plot_dc <- dot_plot_dc + theme(legend.position = "none")


# Extended dot plot DC
Idents(seurat) <- "annotation_paper"
goi_aDC <- c("XCR1", "CLEC9A", "ANPEP", "FBXO27", "CLEC10A",
             "SIRPA", "DENND3", "IL3RA", "JCHAIN", "CLEC4C",
             "LAMP3", "CCR7", "AIRE", "FOXD4", "TNFRSF11A",
             "ST7", "CST7", "TNFRSF11B", "NEK6",
             "MT2A", "CXCL9", "CXCL10", "CCL22", "CCL17",
             "CCL19", "CD40", "CD80", "CD274", "PDCD1LG2",
             "CD74", "HLA-DRA")
aDC_levels <- c("aDC1", "aDC2", "aDC3")
aDC <- seurat[, seurat$annotation_20220215 %in% aDC_levels]
dot_plot_aDC <- DotPlot(aDC, features = rev(goi_aDC)) +
  coord_flip() +
  scale_color_distiller(palette = "Blues", direction = 1) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1))
legend_dot_plot <- as_ggplot(get_legend(dot_plot_aDC))
dot_plot_aDC <- dot_plot_aDC + theme(legend.position = "none")


# Violin plots DC signatures across aDC
markers_dc_list <-  purrr::map(dc_levels, function(x) {
  df <- FindMarkers(dc, ident.1 = x, only.pos = TRUE)
  df <- df %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < 0.001) %>%
    arrange(desc(avg_log2FC))
  markers <- df$gene[1:50]
  markers
})
names(markers_dc_list) <- dc_levels
aDC <- UCell::AddModuleScore_UCell(aDC, features = markers_dc_list)
vars <- c("DC1.precursor_UCell", "DC1.mature_UCell", "DC2_UCell", "DC3_UCell",
          "DC4_UCell", "DC5_UCell", "IL7R.DC_UCell")
aDC$annotation_paper <- droplevels(aDC$annotation_paper)
vlns_aDC_l <- map(vars, function(x) {
  p <- VlnPlot(
    aDC,
    x,
    group.by = "annotation_paper",
    cols = colors_myel[levels(aDC$annotation_paper)],
    pt.size = 0
  ) &
    NoLegend() &
    stat_summary(
      fun = "median",
      geom = "point",
      color = "black",
      size = 0.75
    ) &
    ylab("UCell Score") &
    theme(
      axis.title.x = element_blank(),
      axis.text = element_text(size = 6, color = "black"),
      axis.title.y = element_text(size = 6),
      plot.title = element_text(size = 8, hjust = 0.5)
    )
  p
})
vlns_aDC_arranged <- wrap_plots(vlns_aDC_l, nrow = 3, ncol = 3)

# Dot plot that shows that MMP12, C1QA and SELENOP are
# specific to slancytes across all the tonsil atlas


# Plot proportion of clusters across donors
ids <- str_remove(unique(seurat$donor_id), "BCLL-")
ids <- as.numeric(str_remove(ids, "-T$"))
seurat$donor_id <- factor(
  seurat$donor_id,
  levels = unique(seurat$donor_id)[order(ids)]
)
proportion_df <- seurat@meta.data %>%
  group_by(donor_id, annotation_paper) %>%
  summarise(n_cells = n()) %>%
  mutate(
    n_cells_total = sum(n_cells),
    percentage_cells = round(n_cells / n_cells_total * 100, 3)
  )
proportion_gg <- proportion_df %>%
  ggplot(aes(donor_id, percentage_cells, fill = annotation_paper)) +
    geom_col() +
    scale_fill_manual(values = colors_myel) +
    labs(x = "donor ID", y = "Percentage of cells (%)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 6),
      axis.title = element_text(size = 6),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(override.aes = list(size = 0.1)))
proportion_gg_legend <- as_ggplot(get_legend(proportion_gg))


# Save
dot_plots <- ggarrange(dot_plot_dc, dot_plot_aDC, common.legend = TRUE)
ggsave(
  filename =  here("results/paper/figures/supplementary_figure_DC_aDC_dotplots.pdf"),
  plot = dot_plots,
  device = cairo_pdf,
  width = 18, 
  height = 14.5, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/supplementary_figure_myeloid_violin_plots_aDC.pdf"),
  plot = vlns_aDC_arranged,
  device = cairo_pdf,
  width = 10, 
  height = 9, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/supplementary_figure_myeloid_proportions.pdf"),
  plot = proportion_gg,
  device = cairo_pdf,
  width = 10, 
  height = 10, 
  units = "cm"
)
ggsave(
  filename =  here("results/paper/figures/supplementary_figure_myeloid_proportions_legend.pdf"),
  plot = proportion_gg_legend,
  device = cairo_pdf,
  width = 15, 
  height = 10, 
  units = "cm"
)
