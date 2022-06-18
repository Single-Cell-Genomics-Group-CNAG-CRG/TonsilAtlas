# header
# test to track changes new
scatterplot <- function(x, y, color) {
  df <- data.frame(
    feat1 = se_sub@assays$MAGIC_Spatial@data[x, ],
    feat2 = se_sub@assays$MAGIC_Spatial@data[y, ],
    color = se_sub@meta.data[, color]
  )
  
  ggplot2::ggplot(df,
                  ggplot2::aes(x = feat1,
                               y = feat2,
                               color = color)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = x,
      y = y,
      color = color)
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

convert_symb_entrez <- function(gene_vec,
                                annotation){
  # This function converts a vector of gene symbols to its entrezID form so that it can be used to perform enrichment
  ###
  # gene_vec: vector of gene symbols to conver to entrezID
  ###
  # Returns a vecotr of entrezIDs
  ###
  
  # library(org.Hs.eg.db)
  if (annotation == "org.Hs.eg.db") {
    require(org.Hs.eg.db)
    # Convert symbols to ENTREZID
    DE_entrezid <- AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                         keys = unlist(gene_vec),
                                         column = "ENTREZID",
                                         keytype = "SYMBOL")
    
  } else if (annotation == "org.Mm.eg.db") {
    
    require(org.Mm.eg.db)
    # Convert symbols to ENTREZID
    DE_entrezid <- AnnotationDbi::mapIds(x = org.Mm.eg.db,
                                         keys = unlist(gene_vec),
                                         column = "ENTREZID",
                                         keytype = "SYMBOL")
    
  }
  
  return(DE_entrezid)
  
}


#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

gene_enrichment_GO <- function(gene_de, 
                               gene_universe,
                               pvalueCutoff = 0.05,
                               testDirection = "over",
                               annotation = "org.Hs.eg.db",
                               ontology = "BP") {
  ####
  # This function performs a hyper geometric test with a set of differentially expressed genes over a defined gene univers. 
  # It is set to be used for human genes only.
  # It requires the packages org.Hs.eg.db and GOstats to be installed.
  ####
  # gene_universe: vector with all the genes evaluated in order to set the gene universe. Must be ENTREZID names for the genes.
  # gene_de: vector with all the differentially expressed genes from the gene universe. Must be ENTREZID names for the genes.
  # pvalueCutoff: Cutoff to use for the Pvalue
  # testDirection: A string which can be either "over" or "under". This determines whether the test performed detects over or under represented GO terms.
  # gene_enrichment_GO: which ontology do we want to test BP (biological processes), CC (cellular component), or MF (molecular function)
  ####
  # Returns: hgOver with the enriched GO terms 
  ####
  # library(org.Hs.eg.db) # It is called in convert_symb_entrez
  library(GOstats)
  
  # Convert Target genes from symbol to entrezID
  target_enrich <- convert_symb_entrez(gene_vec = unique(gene_de),
                                       annotation = annotation)
  
  # Convert Universe genes from symbol to entrezID
  univers_g <- convert_symb_entrez(gene_vec = unique(gene_universe),
                                   annotation = annotation) 
  
  # Create a GOHyperGParams instance
  params <- new("GOHyperGParams",
                geneIds = target_enrich, 
                universeGeneIds = univers_g,
                annotation = annotation,
                ontology = ontology,
                pvalueCutoff = pvalueCutoff,
                testDirection = testDirection)
  
  # Carry out hyper geometric test
  hgOver <- GOstats::hyperGTest(params)
  return(hgOver)
  
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

GO_visualization <- function(ds) {
  plt_ls <- lapply(unique(ds$cluster), function(i) {
    tmp_plt <- ds %>%
      dplyr::filter(cluster == i) %>%
      dplyr::arrange(desc(OddsRatio)) %>%
      head(20) %>% 
      ggplot(.) +
      geom_point(aes(x = OddsRatio,
                     y = reorder(Term, OddsRatio),
                     size = -Pvalue,
                     color = Pvalue)) +
      labs(title = sprintf("Cluster %s", i - 1)) +
      scale_color_gradient(low = "green",
                           high = "red") +
      theme_classic()
    return(tmp_plt)
  })
  return(plt_ls)
}

################################################################################
################################################################################
################################################################################

qc_hist_plots <- function(se,
                          nfeat = "nFeature_Spatial",
                          ncount = "nCount_Spatial",
                          slot = "counts",
                          assay = "Spatial",
                          percent.mito = NULL,
                          percent.ribo = NULL) {
  ### se: Seurat object
  require(ggplot2)
  require(Seurat)
  require(ggpubr)
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se[[]], 
                            ggplot2::aes_string(nfeat),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::ggtitle("Unique genes per spot") +
    ggplot2::labs(x = "Number of Detected Genes",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se[[]],
                            ggplot2::aes_string(ncount),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::ggtitle("Total counts per spots") +
    ggplot2::labs(x = "Library Size (total UMI)",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  count_mtrx <- Seurat::GetAssayData(object = se, slot = slot, assay = assay)
  gene_attr <- data.frame(nUMI = Matrix::rowSums(count_mtrx), 
                          nSpots = Matrix::rowSums(count_mtrx > 0))
  p3 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = gene_attr,
                            ggplot2::aes(nUMI),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::scale_x_log10() +
    ggplot2::ggtitle("Total counts per gene (log10 scale)") +
    ggpubr::theme_pubr()
  
  p4 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = gene_attr,
                            ggplot2::aes(nSpots),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::ggtitle("Total spots per gene") +
    ggpubr::theme_pubr()
  
  # Collect all genes coded on the mitochondrial genome
  if (is.null(percent.mito)) {
    mt.genes <- grep(pattern = "^MT-",
                     x = rownames(count_mtrx),
                     value = TRUE,
                     ignore.case = TRUE)
    se_obj[["percent.mito"]] <- (Matrix::colSums(count_mtrx[mt.genes, ]) /
                                   Matrix::colSums(count_mtrx)) * 100
  } else {
    se_obj[["percent.mito"]] <- se_obj[[percent.mito]]
  }
  
  
  # Collect all genes coding for ribosomal proteins
  if (is.null(percent.ribo)) {
    rp.genes <- grep(pattern = "^RPL|^RPS",
                     x = rownames(count_mtrx),
                     value = TRUE,
                     ignore.case = TRUE)
    se_obj[["percent.ribo"]] <- (Matrix::colSums(count_mtrx[rp.genes, ]) /
                                   Matrix::colSums(count_mtrx)) * 100
  } else {
    se_obj[["percent.ribo"]] <- se_obj[[percent.ribo]]
  }
  
  p5 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se_obj[[]], 
                            ggplot2::aes(percent.mito),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::ggtitle("Mitochondrial % per spot") +
    ggplot2::labs(x = "Mitochondrial % ",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  p6 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = se_obj[[]], 
                            ggplot2::aes(percent.ribo),
                            fill = "red",
                            alpha = 0.7,
                            color = "red",
                            bins = 50) +
    ggplot2::ggtitle("Ribosomal % per spot") +
    ggplot2::labs(x = "Ribosomal % ",
                  y = "Number of Spots") +
    ggpubr::theme_pubr()
  
  
  return(list(p1, p2, p3, p4, p5, p6))
}

################################################################################
################################################################################
################################################################################

qc_st_plots <- function(se,
                        nfeat = "nFeature_Spatial",
                        ncount = "nCount_Spatial",
                        slot = "counts",
                        assay = "Spatial",
                        percent.mito = NULL,
                        percent.ribo = NULL) {
  ### se: Seurat object
  require(ggplot2)
  require(Seurat)
  
  p1 <- Seurat::SpatialFeaturePlot(object = se,
                                   features = nfeat) +
    ggplot2::ggtitle("Unique genes per spot")
  
  p2 <- Seurat::SpatialFeaturePlot(object = se,
                                   features = ncount) +
    ggplot2::ggtitle("Total counts per spots")
  
  count_mtrx <- Seurat::GetAssayData(object = se,
                                     slot = slot, assay = assay)
  gene_attr <- data.frame(nUMI = Matrix::rowSums(count_mtrx), 
                          nSpots = Matrix::rowSums(count_mtrx > 0))
  
  # Collect all genes coded on the mitochondrial genome
  if (is.null(percent.mito)) {
    mt.genes <- grep(pattern = "^MT-",
                     x = rownames(count_mtrx),
                     value = TRUE,
                     ignore.case = TRUE)
    se_obj[["percent.mito"]] <- (Matrix::colSums(count_mtrx[mt.genes, ]) /
                                   Matrix::colSums(count_mtrx)) * 100
  } else {
    se_obj[["percent.mito"]] <- se_obj[[percent.mito]]
  }
  
  
  # Collect all genes coding for ribosomal proteins
  if (is.null(percent.ribo)) {
    rp.genes <- grep(pattern = "^RPL|^RPS",
                     x = rownames(count_mtrx),
                     value = TRUE,
                     ignore.case = TRUE)
    se_obj[["percent.ribo"]] <- (Matrix::colSums(count_mtrx[rp.genes, ]) /
                                   Matrix::colSums(count_mtrx)) * 100
  } else {
    se_obj[["percent.ribo"]] <- se_obj[[percent.ribo]]
  }
  
  p3 <- Seurat::SpatialFeaturePlot(object = se,
                                   features = "percent.mito") +
    ggplot2::ggtitle("Mitochondrial % per spot")
  
  p4 <- Seurat::SpatialFeaturePlot(object = se,
                                   features = "percent.ribo") +
    ggplot2::ggtitle("Ribosomal % per spot")
  
  
  return(list(p1, p2, p3, p4))
}

################################################################################
################################################################################
################################################################################

qc_covar_plots <- function(se,
                           nfeat = "nFeature_Spatial",
                           ncount = "nCount_Spatial",
                           slot = "counts",
                           assay = "Spatial",
                           percent.mito = NULL,
                           percent.ribo = NULL) {
  ### se: Seurat object
  require(ggplot2)
  require(Seurat)
  require(ggpubr)
  
  ## Compute mito and ribo %
  count_mtrx <- Seurat::GetAssayData(object = se,
                                     slot = slot,
                                     assay = assay)
  
  # Collect all genes coded on the mitochondrial genome
  if (is.null(percent.mito)) {
    mt.genes <- grep(pattern = "^MT-",
                     x = rownames(count_mtrx),
                     value = TRUE,
                     ignore.case = TRUE)
    se_obj[["percent.mito"]] <- (Matrix::colSums(count_mtrx[mt.genes, ]) /
                                   Matrix::colSums(count_mtrx)) * 100
  } else {
    se_obj[["percent.mito"]] <- se_obj[[percent.mito]]
  }
  
  
  # Collect all genes coding for ribosomal proteins
  if (is.null(percent.ribo)) {
    rp.genes <- grep(pattern = "^RPL|^RPS",
                     x = rownames(count_mtrx),
                     value = TRUE,
                     ignore.case = TRUE)
    se_obj[["percent.ribo"]] <- (Matrix::colSums(count_mtrx[rp.genes, ]) /
                                   Matrix::colSums(count_mtrx)) * 100
  } else {
    se_obj[["percent.ribo"]] <- se_obj[[percent.ribo]]
  }
  
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = se[[]], 
                        ggplot2::aes_string(x = ncount,
                                            y = nfeat,
                                            color = "percent.mito"),
                        alpha = 0.7) +
    scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
    ggplot2::ggtitle("Library size vs Detected genes") +
    ggplot2::labs(x = "Library Size", y = "Unique UMIs") +
    ggpubr::theme_pubr()
  
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = se[[]], 
                        ggplot2::aes_string(x = ncount,
                                            y = nfeat,
                                            color = "percent.ribo"),
                        alpha = 0.7) +
    ggplot2::scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
    ggplot2::ggtitle("Library size vs Detected genes") +
    ggplot2::labs(x = "Library Size", y = "Unique UMIs") +
    ggpubr::theme_pubr()
  
  
  p3 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = se_obj[[]], 
                        ggplot2::aes_string(x = "percent.mito",
                                            y = nfeat,
                                            color = ncount),
                        alpha = 0.7) +
    ggplot2::scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
    ggplot2::ggtitle("Mitochondrial % per spot") +
    ggplot2::labs(x = "Mitochondrial % ", y = "Unique UMIs") +
    ggpubr::theme_pubr()
  
  p4 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = se_obj[[]], 
                        ggplot2::aes(x = percent.mito,
                                     y = percent.ribo),
                        fill = "red",
                        alpha = 0.7,
                        color = "red") +
    ggplot2::ggtitle("Mitochondrial % per spot") +
    ggplot2::labs(x = "Mitochondrial % ", y = "Ribosomal %") +
    ggpubr::theme_pubr()
  
  return(list(p1, p2, p3, p4))
}


################################################################################
################################################################################
################################################################################

GO_function <- function(marker_ls, univ_mouse, univ_human) {
  
  GO_ls <- lapply(names(marker_ls), function(nm) {
    # print(nm)
    markers_mouse <- marker_ls[[nm]][stringr::str_detect(string = marker_ls[[nm]],
                                                         pattern = "mm10---")]
    
    markers_mouse <- stringr::str_remove(string = markers_mouse,
                                         pattern = "mm10---")
    # print(markers_mouse)
    if (length(markers_mouse) > 1) {
      mouse_entrez <- convert_symb_entrez(
        gene_vec = markers_mouse,
        annotation = "org.Mm.eg.db")
    } else {
      mouse_entrez <- NULL
    }
    
    
    if (!is.null(mouse_entrez)) {
      # print("mouse")
      go_clust_mouse <- gene_enrichment_GO(
        gene_de = markers_mouse,
        gene_universe = univ_mouse,
        pvalueCutoff = 0.05,
        testDirection = "over",
        annotation = "org.Mm.eg.db",
        ontology = "BP"
      )
      go_clust_mouse_df <- summary(go_clust_mouse)
      # print("mouse GO end 1")
      go_clust_mouse_df <- go_clust_mouse_df %>% 
        data.frame(.)
      # print("mouse GO end 2")
      go_clust_mouse_df <- go_clust_mouse_df %>% 
        dplyr::mutate(
          cluster = nm,
          organism = "mouse")
      # print("mouse GO end 3")
      # print(head(go_clust_mouse_df))
    }
    
    markers_human <- marker_ls[[nm]][stringr::str_detect(string = marker_ls[[nm]],
                                                         pattern = "GRCh38-")]
    
    markers_human <- stringr::str_remove(string = markers_human,
                                         pattern = "GRCh38-")
    # print(markers_human)
    if ( length(markers_human) > 1) {
      human_entrez <- convert_symb_entrez(
        gene_vec = markers_human,
        annotation = "org.Hs.eg.db")
    } else {
      human_entrez <- NULL
    }
    
    if (!is.null(human_entrez)) {
      # print("human")
      go_clust_human <- gene_enrichment_GO(
        gene_de = markers_human,
        gene_universe = univ_human,
        pvalueCutoff = 0.05,
        testDirection = "over",
        annotation = "org.Hs.eg.db",
        ontology = "BP"
      )
      go_clust_human_df <- summary(go_clust_human)
      # print(go_clust_human)
      # print(head(go_clust_human_df))
      # print(class(go_clust_human))
      # print(class(go_clust_human_df))
      # print("human GO end 1")
      go_clust_human_df <- go_clust_human_df %>% 
        data.frame()
      # print("human GO end 2")
      go_clust_human_df <- go_clust_human_df %>% 
        dplyr::mutate(
          cluster = nm,
          organism = "human")
      # print("human GO end 3")
      # print(head(go_clust_human_df))
    }
    
    
    if (!is.null(human_entrez) &
        !is.null(mouse_entrez)) {
      go_clust <- dplyr::bind_rows(go_clust_mouse_df,
                                   go_clust_human_df)
    } else if (is.null(human_entrez) &
               !is.null(mouse_entrez)) {
      go_clust <- go_clust_mouse_df
    } else if (!is.null(human_entrez) &
               is.null(mouse_entrez)) {
      go_clust <- go_clust_human_df
    }
    
    return(go_clust)
  })
  # print(GO_ls)
  GO_df <- dplyr::bind_rows(GO_ls)
  
  return(GO_df)
}

################################################################################
################################################################################
################################################################################


################################################################################
################################################################################
################################################################################


################################################################################
################################################################################
################################################################################
