# This script contains the utilities needed to plot figure 4 (cytotoxic)


# COLORS
colors_rna <- c(
  "Naive CD8 T" = "#ffb6c1",
  "SCM CD8 T" = "#ffc4a6",
  "CM CD8 T" = "#e46c70",
  "RM CD8 T" = "#ce262f",
  "RM CD8 activated T" = "#862222",
  "CD8 Tf" = "#74b2d0",
  "DC recruiters CD8 T" = "#9E8073", 
  "IFN CD8 T" = "#441815",
  "Nksig CD8 T" = "#E58F47",
  
  "CD56+ gd T" = "#a2c8c7",
  "TCRVδ+ gd T" = "#5ea19e",
  "MAIT" = "#008080",
  
  "CD16-CD56+ NK" = "#a19d41",
  "CD16-CD56- NK" = "#fcf75e",
  "CD16+CD56- NK" = "#cec94f",
  
  "ILC1" = "#80eda8",
  "NKp44+ ILC3" = "#6ac18a",
  "NKp44- ILC3" = "#55976d",
  
  "DN"= "#003b59" 
  
)
goi_rna <- c("CD3D", "CD8A", "TRDC", "TRGC2", "KLRB1", "KLRC1", "KLRD1",
             "KLRF1", "GNLY", "IL7R", "IL1R1", "HAVCR2", "NOSIP",
             "LEF1", "TSHZ2", "FAM117B", "ANXA1", "S100A4", "CD69", "CCL5",
             "GZMK", "CXCR6", "CXCL13", "PDCD1", "ICOS", "TIGIT", "XCL1", "XCL2",
             "IFI44L", "KLRC2", "IL18RAP", "NCAM1", "IKZF3", "FCGR3A", "PRF1",
             "CD200R1", "NCR2", "IL23R", "RORC", "ISG20")

goi_rna_supplement <- list(
  "Naive CD8 T" = c("CCR7", "BACH2"),
  "SCM CD8 T" = c("PLAC8", "ATP8A1", "TCF7"),
  "CM CD8 T" = c("S100A4"),
  "RM CD8 T" = c("EOMES", "ZEB2", "ITGA1", "ITGAE", "TOX"),
  "RM CD8 activated T" = c("GZMH", "CCL5", "GZMK", "GZMA", "CXCR6", "HLA-DRB1", "HLA-DPA1"),
  "CD8 Tf" = c("CD200", "LAG3", "IFNG"),
  "DC recruiters CD8 T" = c("CCL4", "CD99"),
  "IFN CD8 T" = c("IFIT1", "IFI44", "IFIT3", "IFI6"),
  "Nksig CD8 T"= c("NKG7", "CX3CR1"),
  "CD56+ gd T" = c("ZNF683"),
  "TCRVδ+ gd T" = c(),
  "MAIT" = c("RORA"),
  "CD16-CD56+ NK" = c("SELL", "IFNG-AS1"),
  "CD16-CD56- NK" = c("CD160", "CD96"),
  "CD16+CD56- NK" = c("KLF2", "FGFBP2", "CX3CR1", "TBX21"),
  "NKp44+ ILC3" =  c("IL4I1", "AHR"),
  "DN" = c("HAVCR2", "CD38", "GPR183")
)

