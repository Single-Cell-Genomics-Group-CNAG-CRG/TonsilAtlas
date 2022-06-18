# This script contains parameters, variables and functions used to create
# figure 3


# FUNCTIONS
theme_nothing2 <- function() {
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
}
change_dot_size <- function(size) {
  guides(color = guide_legend(override.aes = list(size = size)))
}


# PATHS
path_to_epi <- "~/Google Drive/My Drive/single_cell/PhD/B_cell_atlas/tonsil_atlas/tonsil_atlas_annotation/current/level_4/epithelial/epithelial_annotated_level_4.rds"
path_to_myel <- "~/Google Drive/My Drive/single_cell/PhD/B_cell_atlas/tonsil_atlas/tonsil_atlas_annotation/current/level_5/myeloid/myeloid_annotated_level_5.rds"
path_to_spatial_epith <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/spatial_transcriptomics/results/R_objects/tonsil_integrated_annotated.rds"
path_to_save <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/plots/paper/non_lymphoid/epithelial_umap_and_dot_plot.pdf"


# COLORS
colors_epi <- c(
  "Basal cells" = "#1b9e77",
  "VEGFA+" = "#d95f02",
  "Surface epithelium" = "#7570b3",
  "Outer surface" = "#e7298a",
  "Crypt" = "#66a61e",
  "FDCSP epithelium" = "#e6ab02"
)
colors_myel <- c(
  "DC1 precursor" = "#ff99d6",
  "DC1 mature" = "#bf408c",
  "DC2" = "#f99fa7",
  "DC3" = "#e06c6c",
  "DC4" = "#8f523d",
  "DC5" = "#d73919",
  "IL7R DC" = "#e1b40e",
  "aDC1" = "#e2f155",
  "aDC2" = "#b4d84b",
  "aDC3" = "#6ead34",
  "M1 Macrophages" = "#3f56ca",
  "Monocytes" = "#cd87de",
  "Mast cells" = "#825e7b",
  "Neutrophil Granulocytes" = "#aa7c60",
  "Cycling" = "gray",
  "IL7R MMP12 macrophages" = "#91f8ca",
  "C1Q HLA macrophages" = "#68cde8",
  "SELENOP FUCA1 PTGDS macrophages" = "#519bdb",
  "ITGAX ZEB2 macrophages" = "#7ef1eb"
)

# GENES OF INTEREST (GOI)
goi_epi <- list(
  Basal_cells = c("KRT5", "KRT14", "S100A2"),
  VEGFA = c("VEGFA", "MIR205HG", "PAX1"),
  Surface_epithelium = c("KRT4", "KRT13", "KRT78", "KRT80", "MAL", "SPRR3", "TMPRSS2", "TMPRSS11B"),
  Outer_surface = c("LCE3A", "LCE3D", "LCE3E", "SPRR2D", "SPRR2E", "CNFN"),
  Crypt = c("KRT8", "CD63", "IL1B", "IFI27", "S100A6", "SPIB", "MARCKSL1"),
  FDCSP_epithelium = c("FDCSP", "KRTDAP", "CALML5")
)
goi_myeloid <- c("CLEC9A", "XCR1", "CD1C", "VCAN", "FCGR3A", "SIGLEC6",
                 "IL7R", "CCR7", "AIRE", "CCL19", "TNIP3", "CXCL8", "S100A12",
                 "TPSAB1", "TOP2A")
# goi_myel <- list(
#   "DC1 precursor" = c("CLEC9A", "CADM1"),
#   "DC1 mature" = c("XCR1"),
#   "DC2" = c("CLEC10A", "CD1C", "HLA-DQB1", "HLA-DPB1"),
#   "DC3" = c("VCAN", "S100A8", "S100A9"),
#   "DC4" = c("SERPINA1", "FTL", "FCGR3A"),
#   "DC5" = c("SIGLEC6", "AXL","CD2", "TCF4", "LILRA4"),
#   "IL7R DC" = c("IL7R"),
#   "aDC1" = c(),
#   "aDC2" = c(),
#   "aDC3" = c(),
#   "IL7R MMP12 macrophages" = c("SLAMF7", "MMP14", "MMP25", "MMP12", "MMP9"),
#   "ITGAX ZEB2 macrophages" = c("SCARB1", "SCARB2", "LY75"),
#   "C1Q HLA macrophages" = c("C1QA", "C1QB", "HLA-DRA"),
#   "SELENOP FUCA1 PTGDS macrophages" = c("APOC1", "APOE", "FUCA1"),
#   "M1 Macrophages" = c(),
#   "Monocytes" = c("S100A4"),
#   "Mast cells" = c("TPSAB1", "TPSB2"),
#   "Granulocytes" = c("CXCL8", "PI3"),
#   "Cycling" = c()
# )
# goi_myel_spatial <- c("MMP14", "C1QA", "SELENOP", "SCARB1")

