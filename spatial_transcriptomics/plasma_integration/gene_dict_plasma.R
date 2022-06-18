# plasma_markers <- "{plasma}/{robj_dir}/Markers_interest_PC.xlsx" %>%
#   glue::glue() %>%
#   here::here() %>%
#   xlsx::read.xlsx(
#     file = .,
#     sheetName = "Hoja1")
# 
# signatures <- c("LZ", "LZ.PRDM1.PC", "LZ.PRDM1.XBP1.IRF4.PC", "DZ", "DZ.precPC",
#   "csMBC", "MBC.early.PC.prec", "MBC..PC.prec")
# 
# plasma_genes <- plasma_markers %>%
#   dplyr::select(-Comment) %>% 
#   tidyr::pivot_longer(
#     cols = signatures,
#     names_to = "signature",
#     values_to = "expression") %>%
#   dplyr::filter(expression > 0) %>%
#   dplyr::group_split(signature) %>%
#   purrr::map(~ dplyr::pull(., "Gene"))
#   
# names(plasma_genes) <- sort(unique(signatures))
# plasma_genes[["DZ_2"]] <- c("OAZ1", "AICDA", "H3", "MKI67", "POLH", "TOP2A")
# plasma_genes[["LZ_2"]] <- c("LAG3", "ITGB8", "PDCD1", "TIGIT", "BCL2", "PECAM1",
#   "LY6E", "B7-H3", "HLA-DRB1", "PSMB10", "TNF", "ARG1", "HLA-E", "STAT1")
# 
# plasma_genes[["IGH"]] <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE", "IGHM")
# plasma_genes[["interest"]] <- c("XBP1", "CREB3L2", "CD9", "CD44", "PRDM1", "IRF4")

# - **Dark Zone**: 
# - **Light Zone**: LAG3, ITGB8, PDCD1, TIGIT, BCL2, PECAM1, LY6E, B7-H3 (CD276), HLA-DRB1, PSMB10, TNF, ARG1, HLA-E, STAT1

## Genes of interest
plasma_genes <- list()
plasma_genes[["LZ"]] <- c("SUGCT", "AICDA","CXCR4")
plasma_genes[["DZ"]] <- c("LMO2", "CD83", "BCL2A1")
plasma_genes[["GCBC"]] <- c("BCL6", "IRF8", "MEF2B", "MS4A1", "PAX5")
plasma_genes[["PC"]] <- c("PRDM1", "XBP1", "IRF4", "SLAMF7", "SSR4", "MZB1", "DERL3", "CREB3L2")
plasma_genes[["IGHG"]] <- c("IGHG1","IGHG2", "IGHG3","IGHG4", "IGHA1", "IGHA2","IGHM","IGHD", "IGHV3-20","IGHV3-43")
plasma_genes[["(Im)Mature"]] <- c("CD9", "CD44")
plasma_genes[["Memory"]] <- c("BANK1", "CELF2", "TXNIP")
# plasma_genes[["Proliferation"]] <- c("H2AFZ", "HMGN2", "MCM6", "MKI67", "PCMA", "SLAMF7", "TOP2A", "TUBA1B")
plasma_genes[["Proliferation"]] <- c("MKI67", "TOP2A", "TUBA1B")
plasma_vec <- unique(unlist(plasma_genes))
