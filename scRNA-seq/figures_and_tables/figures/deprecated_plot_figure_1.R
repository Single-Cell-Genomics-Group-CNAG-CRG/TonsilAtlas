# This script plots the figure 1


# Load packages
library(ComplexHeatmap)
library(tidyverse)
library(ggpubr)
library(ggmap)
library(ggrastr)
library(patchwork)


# Source utils
path_to_utils <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/bin/utils_figure1.R"
source(path_to_utils)


# Read data
figure_1_df <- read_csv(path_to_umap_df, col_types = "cccddccccd")
donor_metadata <- read_csv(path_to_donor_metadata)
spatial_transcript_metadata <- read_csv(path_to_spatial_metadata)


# Plot
figure_1_df_sub <- figure_1_df %>%
  dplyr::filter(annotation_prob > 0.6 | is.na(annotation_prob))
umap_title <- format(nrow(figure_1_df_sub), big.mark = ",", scientific = FALSE)
umap_all <- figure_1_df_sub %>%
  ggplot(aes(UMAP_1_level_1, UMAP_2_level_1, color = annotation_figure_1)) +
    geom_point(size = 0.25) +
    theme_classic() +
    scale_color_manual(values = colors_figure_1) +
    labs(title = str_c(umap_title, "cells", sep = " "), color = "") +
    theme_nothing2() +
    theme(
      plot.title = element_text(size = 8, hjust = 0),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.text.align = 0
    ) +
    guides(color = guide_legend(nrow = 4, byrow = TRUE, override.aes = list(size = 2))) +
    coord_fixed(ratio = 1)
legend_all <- as_ggplot(get_legend(umap_all))
umap_all <- umap_all + no_legend()
umap_all <- rasterize(umap_all, dpi = 300)


# UMAP split by assay
assays_df$n_cells <- purrr::map_int(assays_df$assay, function(x) {
  out <- sum(figure_1_df_sub$assay == x)
  out
})
umaps_list <- purrr::map(assays_df$assay, function(x) {
  df <- figure_1_df_sub[figure_1_df_sub$assay == x, ]
  n_cells <- format(
    assays_df[assays_df$assay == x, "n_cells"],
    big.mark = ",",
    scientific = FALSE
  )
  p <- df %>%
    ggplot(aes(UMAP_1_level_1, UMAP_2_level_1)) +
      geom_point(size = 0.1, color = assays_df[assays_df$assay == x, "color"]) +
      labs(title = x, subtitle = str_c(n_cells, "cells", sep = " ")) +
      theme_classic() +
      labs(color = "") +
      theme_nothing2() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 7.5),
        plot.subtitle = element_text(size = 6.5),
        legend.position = "bottom",
        legend.background = element_blank()
      ) +
      guides(color = guide_legend(override.aes = list(size = 2))) +
      coord_fixed(ratio = 1)
  rasterize(p, dpi = 300)
})
umaps_assay_arr <- (umaps_list[[1]] | umaps_list[[2]]) /
                   (umaps_list[[3]] | umaps_list[[4]])


# Heatmap to show donor metadata
table_assays <- table(figure_1_df$assay, figure_1_df$donor_id)
table_assays <- unclass(table_assays > 0) * 1
donor_metadata <- as.data.frame(donor_metadata)
rownames(donor_metadata) <- donor_metadata$donor_id
donor_metadata <- donor_metadata[, c("age_group", "sex", "hospital")]
table_assays <- table_assays[assays, rownames(donor_metadata)]
spatial_trans_row <- ifelse(
  colnames(table_assays) %in% spatial_transcript_metadata$donor_id,
  1,
  0
)
table_assays <- rbind(table_assays, spatial_trans_row)
rownames(table_assays)[rownames(table_assays) == "spatial_trans_row"] <- "Spatial"
# Note: in a CITE-seq experiment, we multiplexed donors BCLL-2-T, BCLL-6-T, BCLL-10-T y BCLL-12-T
table_assays["CITE-seq", c("BCLL-2-T", "BCLL-6-T", "BCLL-10-T", "BCLL-12-T")] <- 2
donor_metadata <- donor_metadata[new_order, ]
table_assays <- table_assays[, new_order]
donor_annotation_heatmap <- HeatmapAnnotation(
  df = donor_metadata,
  col = colors_metadata,
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    age_group = list(
      title = "age group",
      at = c("kid", "young_adult", "old_adult"),
      labels = c("kid", "young adult", "old adult")
    ),
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6),
    direction = "vertical",
    ncol = 1,
    grid_height = unit(0.25, "cm"),
    grid_width = unit(0.25, "cm")
  )
)
samples_heatmap <- Heatmap(
  table_assays,
  col = colors_heat,
  rect_gp = gpar(col = "white", lwd = 2),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = donor_annotation_heatmap, 
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title = "Sample",
    at = c(0, 1, 2), 
    labels = c("no", "yes", "multiplexed"),
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6),
    direction = "vertical",
    ncol = 1,
    grid_height = unit(0.25, "cm"),
    grid_width = unit(0.25, "cm")
  )
)
samples_heatmap <- draw(
  samples_heatmap,
  annotation_legend_side = "right",
  heatmap_legend_side = "right"
) %>%
  grid.grabExpr() %>%
  ggplotify::as.ggplot()


# Arrange
top <- (empty_plot | samples_heatmap) +
  patchwork::plot_layout(widths = c(0.35, 0.65))
# mid <- (umap_all | umaps_assay_arr)/ legend_all
# mid <- mid + patchwork::plot_layout(heights = c(0.95, 0.05))
mid <- (umap_all | umaps_assay_arr)
fig1 <- (top / mid) + patchwork::plot_layout(heights = c(0.35, 0.65))


# Save
ggsave(
  plot = fig1,
  filename = path_to_save,
  height = 21,
  width = 20,
  units = "cm"
)
ggsave(
  plot = legend_all,
  filename = path_to_save_legend,
  height = 5,
  width = 10,
  units = "cm"
)


# Post-processing with Inkscape
## 1. Open inkscape, open figure 1 (cntr+O), change document size to A4 (ctrl + shift + D)
## 2. Eliminate all unnecessary white backgrounds (right click + ungroup + delete)
## 3. Align all legends panel b (heatmap) vertically, move whole heatmap to the top-right corner
## 4. Reduce space between heatmap and UMAPs
## 5. Add panel labels (a-f): sans serif, bold, 11. Use guides and ctrl + D (duplicate)
## 6. Add arrows for UMAP axis: shift + F6 --> open fill and stroke dialong (shift + ctrl + F) +
##    stroke style tab + select triangular marker (3 dialog box in Markers:) pointing right. Width =0,165
##    Duplicate arrow (ctrl + D) --> flip 90º (top left, under Object) --> crisp --> group (ctrl + G)
## 7. Import (File --> Import or ctrl + I) legend panel c (UMAP all). Reduce space between: rows, columns,
##    and point/text. Shorten names (finish with .). Place top right of panel c.
