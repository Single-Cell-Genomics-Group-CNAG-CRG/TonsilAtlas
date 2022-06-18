# This script plots the supplementary figure that validates the integration:
# 1. Preservation of biological variability: UMAP King et al colored by their annotation
# 2. Removal technical/unwanted variability: UMAP with all assays (multiome, King, 3P) 
# colored and split by assay. Include before and after integration +
# LISI score


# Load package
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrastr)
library(ggpubr)
library(pals)
library(here)


# Source utils
source(here("scRNA-seq/bin/utils_figure1.R"))


# Define paths
path_to_data <- here("data/raw_data_figures/umaps_rna_batch_correction.csv")
path_to_lisi_scores <- here("data/raw_data_figures/lisi_scores_rna.csv")
path_to_save <- here("results/paper/figures/supplementary_figure_integration.png")
path_to_save_legend <- here("results/paper/figures/supplementary_figure_integration_legend.pdf") 

# Read data
umaps_df <- read_delim(path_to_data, delim = ";", col_names = TRUE)
supp_lisi_df <- read_delim(path_to_lisi_scores, delim = ";", col_names = TRUE)


# Plot
cols <- cols_assays_df$color[cols_assays_df$assay %in% c("scRNA-seq", "Multiome")]
cols <- c(cols, "gray50")
names(cols) <- c("scRNA-seq", "Multiome", "King et al.")
umaps_df$assay <- case_when(
  umaps_df$assay == "3P" ~ "scRNA-seq",
  umaps_df$assay == "5P" ~ "King et al.",
  umaps_df$assay == "multiome" ~ "Multiome"
)
umaps_df$assay <- factor(umaps_df$assay, levels = names(cols))
umaps_df$is_integrated <- factor(
  umaps_df$is_integrated,
  levels = c("unintegrated", "integrated")
)
umaps_batch <- umaps_df %>%
  ggplot(aes(UMAP_1, UMAP_2, color = assay)) +
    geom_point(shape = ".", alpha = 0.6) +
    scale_color_manual(values = cols, breaks = names(cols)) +
    facet_grid(assay ~ is_integrated) +
    labs(x = "UMAP1", y = "UMAP2") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 7, color = "black"),
      axis.title = element_text(size = 6),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(size = 0.25)
    )
# umaps_batch <- rasterize(umaps_batch, dpi = 300)
n_cell_types <- length(unique(umaps_df$annotation_king))
new_levels <- c("Naive", "MBC", "FCRL2/3high GC", "MBC FCRL4+", "Activated", "preGC", "LZ GC", "GC", "DZ GC",
  "Cycling", "prePB", "Plasmablast", "Cycling T", "CD4+", "CD4+ NCM", "TfH", "Treg", "TfR",
  "CD8+ Cytotoxic", "CD8+ NCM", "TIM3+ DN", "NK", "ILC", "cDC1", "MAC1", "MAC2", "MAC3",
  "FDC", "pDC", "Precursors")
umaps_king <- umaps_df %>%
  filter(assay == "King et al." & is_integrated == "integrated")
umaps_king$annotation_king <- factor(umaps_king$annotation_king, new_levels)
umap_preservation <- umaps_king %>%
  ggplot(aes(UMAP_1, UMAP_2, color = annotation_king)) +
    geom_point(size = 0.01) +
    scale_color_manual(values = pals::glasbey(n = n_cell_types)) +
    ggtitle("King et al.") +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.spacing = unit(0, "cm"),
      legend.box.spacing = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.key.height = unit(0, "cm"),
      legend.key.width = unit(0, "cm"),
      plot.title = element_text(hjust = 0.5, size = 7),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.text = element_text(size = 6)
    ) +
    guides(color = guide_legend(override.aes = list(size = 1)))
# umap_preservation <- rasterize(umap_preservation, dpi = 300)


# LISI scores
lisi_vars <- c("lisi_sex", "lisi_age_group", "lisi_is_hashed",
               "lisi_hospital", "lisi_assay")
confounders_vars <- c("sex", "age_group","is_hashed", "hospital", "assay")
cols_integration <- c(unintegrated = "gray60", integrated ="limegreen")
supp_lisi_df$is_integrated <- factor(
  supp_lisi_df$is_integrated,
  levels = names(cols_integration)
)
supp_lisi_df$age_group <- factor(
  supp_lisi_df$age_group,
  levels = c("kid", "young_adult", "old_adult")
)
levels(supp_lisi_df$age_group) <- c("kid", "young adult", "old adult")
supp_lisi_df$is_hashed[supp_lisi_df$assay == "multiome"] <- "not_hashed"
supp_lisi_df$is_hashed[supp_lisi_df$is_hashed == "not_hashed"] <- "not hashed"
supp_lisi_df$assay <- case_when(
  supp_lisi_df$assay == "3P" ~ "scRNA-seq",
  supp_lisi_df$assay == "multiome" ~ "Multiome"
)
supp_lisi_df$assay <- factor(
  supp_lisi_df$assay,
  levels = c("scRNA-seq", "Multiome")
)
supp_lisi_df$hospital <- factor(supp_lisi_df$hospital, levels = c("CIMA", "Clinic"))
levels(supp_lisi_df$hospital) <- c("Pamplona", "Barcelona")
boxplots_lisi <- purrr::map2(confounders_vars, lisi_vars, function(x, y) {
  p <- ggplot(supp_lisi_df, aes_string(x, y, fill = "is_integrated")) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(
      values = cols_integration,
      breaks = names(cols_integration)
    ) +
    labs(title = x, y = "LISI") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 7),
      axis.text = element_text(size = 6, color = "black"),
      axis.title.y = element_text(size = 6), 
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      axis.line = element_line(size = 0.25)
    )
  p
})
boxplots_lisi[[2]] <- boxplots_lisi[[2]] +
  ggtitle("age group") +
  scale_x_discrete(labels = c("kid", "young\nadult", "old\nadult"))
  # theme(axis.text.x = element_blank())
legend_boxplot <- as_ggplot(get_legend(boxplots_lisi[[1]]))
boxplots_lisi <- purrr::map(boxplots_lisi, function(p) {
  p + theme(legend.position = "none")
})
boxplots_arranged <- ggarrange(
  plotlist = boxplots_lisi,
  nrow = 3,
  ncol = 2
)


# Arrange figure
left <- umaps_batch
right <- (boxplots_arranged / umap_preservation) +
  plot_layout(heights = c(3, 2))
fig <- left | right +
  plot_layout(widths = c(2, 3))


# Save
ggsave(
  plot = fig,
  filename = path_to_save,
  width = 20,
  height = 18,
  dpi = 300,
  units = "cm"
)
ggsave(
  plot = legend_boxplot,
  filename = path_to_save_legend,
  width = 20,
  height = 18,
  units = "cm"
)


