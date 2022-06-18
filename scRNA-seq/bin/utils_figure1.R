# This script contains parameters, variables and functions used to create
# figure 1


# FUNCTIONS
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
path_to_umap_df <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/tables/data_figures_paper/umaps_figure_1.csv"
path_to_donor_metadata <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/data/tonsil_atlas_donor_metadata.csv"
path_to_spatial_metadata <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/spatial_transcriptomics/01-spaceranger/data/sample_id.txt"
path_to_save <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/plots/paper/figure1_umaps.pdf"
path_to_save_legend <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/plots/paper/legend_figure1_umap.pdf"


# COLORS
# colors_figure_1 <- c(
#   "NBC_MBC" = "#4075b5",
#   "GCBC" = "#6eb8e2",
#   "PC" = "#9ae9f9",
#   "CD4 T" = "#3eac67",
#   "cycling T cells" = "#59f7a8",
#   "naive CD8 T" = "#e43fde",
#   "CD8 T" = "#a835ac",
#   "TIM3 DN" = "#401c59",
#   "NK" = "#ec6ab6",
#   "ILC" = "#8d4466",
#   "DC" = "#e89064",
#   "aDC" = "#573232",
#   "Macrophages" = "#cb664d",
#   "Monocytes" = "#ff5254",
#   "Mast" = "black",
#   "Granulocytes" = "#c6a23f",
#   "FDC" = "#c9e160",
#   "cycling FDC" = "#7e9a3c",
#   "epithelial" = "#eaef48",
#   "PDC" = "gray30",
#   "pre B/T" = "#3d3b81",
#   "Cycling" = "gray88"
# )

cols_fig1 <- c(
  "NBC" = "#dcf0f4",
  "Activated NBC" = "#7bc6d6", 
  "GCBC" = "#398D9F",
  "MBC" = "#025566", 
  "PC" = "#032872", 
  
  "CD4 T" = "#93331b", 
  "Naive CD4 T" = "#be8476", 
  "Naive CD8 T" = "#e5b7c6", 
  "CD8 T" = "#cb708e", 
  "NK" = "#67253a", 
  "DN" = "#d2a027", 
  "ILC" = "#88657f",
  "cycling T"  = "#cc0000",
  
  "preB/T" = "#4363d8", 
  
  
  "DC" = "#b6dec3", 
  "Mono/Macro" = "#6ac087", 
  "Mast" = "#3a8b55", 
  "Granulocytes" = "#32d732", 
  "cycling myeloid" = "#2b6840",
  
  "FDC" = "#d5ff00", 
  "cycling FDC" = "#bfe600", 
  
  "epithelial" = "#ffdc9f", 
  
  "PDC" = "#f032e6"
)
assays <- c("scRNA-seq", "scATAC-seq", "Multiome", "CITE-seq")
cols_assays <- c("#d52f54", "#73bd56", "#4bafdd", "#f0b635")
cols_assays_df <- data.frame(
  assay = c("scRNA-seq", "scATAC-seq", "Multiome", "CITE-seq"),
  color = c("#d52f54", "#73bd56", "#4bafdd", "#f0b635")
)
colors_metadata <- list(
  sex = c(male = "#00c4aa", female = "#8700f9"),
  age_group = c(kid = "#d08898", young_adult = "#995c96", old_adult = "#563454"),
  hospital = c(Clinic = "#56719f", CIMA = "#c68560")
)
colors_heat <- c("white", "gray70", "gray30")


# MISC
# empty_plot <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#   geom_blank() +
#   theme_nothing()
new_order <- c("BCLL-11-T", "BCLL-12-T", "BCLL-13-T", "BCLL-8-T", "BCLL-9-T",
               "BCLL-10-T", "BCLL-14-T", "BCLL-15-T", "BCLL-6-T", "BCLL-2-T")