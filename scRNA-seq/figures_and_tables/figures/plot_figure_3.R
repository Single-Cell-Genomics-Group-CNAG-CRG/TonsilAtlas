# This script plots figure 3


# Load packages
library(tidyverse)
library(patchwork)
library(ggpubr)
library(Seurat)
set.seed(123)


# Source utils
path_to_utils <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/bin/utils_figure3.R"
source(path_to_utils)


# Read data
epi <- readRDS(path_to_epi)
myel <- readRDS(path_to_myel)
# spatial <- readRDS(path_to_spatial_epith)


###############################################################################
########################### EPITHELIAL (epi) ##################################
###############################################################################

# UMAP
umap_df_epi <- Embeddings(epi, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = factor(
    epi$annotation_level_3,
    levels = names(colors_epi)
  ))
umap_gg_epi <- ggplot(umap_df_epi, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point(size = 1, alpha = 0.85) +
  scale_color_manual(values = colors_epi) +
  theme_classic() +
  theme_nothing2() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85, 0.75),
    legend.text = element_text(size = 6)
  ) +
  change_dot_size(2)


# Dot Plot
goi_epi <- unlist(goi_epi)
names(goi_epi) <- NULL
epi$annotation_figure_3 <- factor(
  epi$annotation_level_3,
  levels = names(colors_epi)
)
dot_plot_epi <- DotPlot(
  epi,
  features = rev(goi_epi),
  group.by = "annotation_figure_3"
)
dot_plot_epi <- dot_plot_epi +
  labs(x = "", y = "") +
  scale_size_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    range = c(0.1, 4.5)
  ) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
dot_plot_epi$guides$colour$title <- "Average\nExpression"
dot_plot_epi$guides$size$title <- "Percent\nExpressed"


###############################################################################
########################### MYELOID (myel) ####################################
###############################################################################

# UMAP
myel$annotation_level_4[myel$annotation_level_4 == "Neutrophil Granulocytes"] <- "Granulocytes" 
myel$annotation_figure_3 <- factor(
  myel$annotation_level_4,
  levels = names(colors_myel)
)
umap_df_myel <- Embeddings(myel, "umap") %>%
  as.data.frame() %>%
  mutate(annotation = myel$annotation_figure_3)
umap_gg_myel <- ggplot(umap_df_myel, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point(size = 1, alpha = 0.85) +
  scale_color_manual(values = colors_myel) +
  theme_classic() +
  theme_nothing2() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 6)
  ) +
  change_dot_size(2)


# Dot Plot
goi_myel <- unlist(goi_myel)
names(goi_myel) <- NULL
dot_plot_myel <- DotPlot(
  myel,
  features = rev(goi_myel),
  group.by = "annotation_figure_3"
)
dot_plot_myel <- dot_plot_myel +
  labs(x = "", y = "") +
  scale_size_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    range = c(0.1, 4.5)
  ) +
  scale_color_distiller(palette = "Blues", direction = 1) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
dot_plot_myel$guides$colour$title <- "Average\nExpression"
dot_plot_myel$guides$size$title <- "Percent\nExpressed"


###############################################################################
############################### ARRANGE #######################################
###############################################################################

fig <- umap | dot_plot
fig  <- fig +
  plot_layout(widths = c(0.67, 0.33))


###############################################################################
############################### SAVE ##########################################
###############################################################################

# Save
ggsave(
  filename = path_to_save,
  plot = fig,
  device = cairo_pdf,
  width = 16, 
  height = 12, 
  units = "cm"
)

