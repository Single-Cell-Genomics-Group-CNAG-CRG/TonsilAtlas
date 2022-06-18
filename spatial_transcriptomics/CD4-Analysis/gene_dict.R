# Dictionary with the genes identifying the different cell types

gene_dict <- list()
gene_dict[["DC"]] <- c("CCR7", "HLA-DRA", "CLEC9A")
gene_dict[["B-cognate"]] <- c("CD40", "CD80", "CD86", "BCL6", "ICOSLG")
gene_dict[["CD4-naive"]] <- c("CCR7", "CD4", "SELPLG", "TCF7", "LEF1", "STAT1", "STAT3")
gene_dict[["CD4-Tfh"]] <- c("BCL6", "CD40LG", "CXCR5", "IL6R", "ICOS", "IL21", "CD28", "SELPLG", "GPR183")
gene_dict[["CD4-non-Tfh"]] <- c("PRDM1", "IL2RA")
gene_dict[["CD4-GC-Tfh"]] <- c("SH2D1A", "CD200", "BTLA", "PDCD1", "CXCR5", "GPR183", "IL4")
gene_dict[["Chemokines"]] <- c("IL2", "IL6", "IL12A", "IL21", "IL10")
gene_dict[["FDC"]] <- c("CR1", "CR2", "CXCL13")
gene_dict[["Proliferation"]] <- c("MKI67", "CDK1", "TOP2A")
gene_dict[["Lymph Endothelial"]] <- c("PECAM1", "CLDN5", "PROX1", "LYVE1")
gene_dict[["T cells"]] <- c("CD3D", "CD3E", "CD8A", "CD8B")
gene_dict[["B cells"]] <- c("MS4A1", "CD79A", "CD40", "ICAM1", "ICOSLG", "CD80")
gene_dict[["DZ"]] <- c("OAZ1", "AICDA", "H3", "MKI67", "POLH")
gene_dict[["LZ"]] <- c("LAG3", "ITGB8", "PDCD1", "TIGIT", "BCL2", "PECAM1",
                       "LY6E", "CD276", "HLA-DRB1", "PSMB10", "TNF", "ARG1",
                       "HLA-E", "STAT1")
gene_dict[["naDC"]] <- c("ITGAX", "XCR1")

gene_vec <- unique(as.character(unlist(gene_dict)))
