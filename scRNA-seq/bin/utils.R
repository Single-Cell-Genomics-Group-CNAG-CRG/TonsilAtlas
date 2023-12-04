# Useful functions used in different scripts


horizontal_barplot <- function(df, categorical_var, continuous_var, ylab) {
  levels_var <- df[[categorical_var]][order(df[[continuous_var]])]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = levels_var)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var)) +
      geom_col() +
      labs(x = "", y = ylab) +
      theme_bw() +
      theme(axis.text = element_text(size = 11),
            axis.title = element_text(size = 13)) +
      coord_flip()
}



customized_boxplot <- function(p) {
  p +
    geom_boxplot(width = 0.75) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 13))
}



horizontal_boxplot <- function(df,
                               categorical_var,
                               continuous_var,
                               fill,
                               ylab,
                               decreasing = FALSE) {
  unordered_lev <- unique(df[[categorical_var]])
  means_cont <- purrr::map_dbl(unordered_lev, function(x) {
    mean(df[[continuous_var]][df[[categorical_var]] == x])
  })
  names(means_cont) <- unordered_lev
  ordered_lev <- names(means_cont)[order(means_cont, decreasing = decreasing)]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = ordered_lev)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var, fill = fill)) +
    geom_boxplot() +
    labs(x = "", y = ylab) +
    theme_bw() +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 13)) +
    coord_flip()
}



horizontal_violin_plot <- function(df,
                               categorical_var,
                               continuous_var,
                               fill,
                               ylab,
                               decreasing = FALSE) {
  unordered_lev <- unique(df[[categorical_var]])
  means_cont <- purrr::map_dbl(unordered_lev, function(x) {
    mean(df[[continuous_var]][df[[categorical_var]] == x])
  })
  names(means_cont) <- unordered_lev
  ordered_lev <- names(means_cont)[order(means_cont, decreasing = decreasing)]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = ordered_lev)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var, fill = fill)) +
    geom_violin() +
    labs(x = "", y = ylab) +
    theme_bw() +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 13)) +
    coord_flip()
}



plot_histogram_qc <- function(df, x, x_lab) {
  df %>%
    ggplot(aes_string(x)) +
    geom_histogram(bins = 100) +
    labs(x = x_lab, y = "Number of Cells") +
    theme_pubr()
}



plot_histogram_doublets <- function(df, x, x_lab, bins) {
  df %>%
    ggplot(aes_string(x)) +
      geom_histogram(bins = bins) +
      xlab(x_lab) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11)
      )
}



plot_density_doublets <- function(df, x, x_lab, color, color_lab) {
  df %>%
    ggplot(aes_string(x = x, color = color)) +
      geom_density() +
      scale_color_brewer(palette = "Dark2") +
      labs(x = x_lab, color = color_lab) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11)
      )
}



plot_boxplot_doublets <- function(df, x, y, fill, y_lab) {
  df %>%
    ggplot(aes_string(x, y, fill = fill)) +
      geom_boxplot() +
      labs(x = "", y = y_lab) +
      scale_fill_brewer(palette = "Dark2") +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 11)
      )
}



plot_scatter_doublets <- function(df, x, y, x_lab, y_lab) {
  df %>%
    ggplot(aes_string(x, y)) +
    geom_point(size = 0.1, alpha = 0.25) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = x_lab, y = y_lab) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}



feature_plot_doublets <- function(seurat_obj, feature) {
  p <- FeaturePlot(
    seurat_obj,
    features = feature,
    cols = viridis::inferno(10),
    pt.size = 0.3
  ) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.text = element_text(size = 11)
    )
  p
}



cluster_diff_resolutions <-  function(seurat_obj,
                                      min_resolution,
                                      max_resolution,
                                      step) {
  # Define resolutions
  resolutions <- seq(from = min_resolution, to = max_resolution, by = step)


  # Find clusters
  resolutions_dfs <- lapply(resolutions, function(x) {
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = x)
    df <- data.frame(
      barcode = colnames(seurat_obj),
      cluster_id = as.character(seurat_obj$seurat_clusters)
    )
    df
  })


  # Merge data frames
  names(resolutions_dfs) <- resolutions
  resolutions_df <- dplyr::bind_rows(resolutions_dfs, .id = "resolution")
  resolutions_df$resolution <- as.numeric(resolutions_df$resolution)
  resolutions_df
}


calculate_oob_accuracy <- function(rf_model) {
  table_rf <- table(rf_model$y, rf_model$predicted)
  oob_acc <- sum(diag(table_rf)) / sum(table_rf) * 100
  oob_acc
}



run_random_forest <- function(resolutions_df,
                              resolution,
                              dim_red_coords,
                              n_cells_per_class = 1000,
                              seed,
                              n_trees,
                              n_cpu,
                              return  = "oob") {

  # Subset to resolution of interest
  resolutions_df <- resolutions_df[resolutions_df$resolution == resolution, ]


  # Get input dataframe and downsample to avoid class imbalance
  input_df <- as.data.frame(dim_red_coords)
  if (all(resolutions_df$barcode == rownames(input_df))) {
    input_df$cluster_id <- factor(resolutions_df$cluster_id)
  } else {
    stop("cell barcodes are different")
  }
  set.seed(seed)
  selected_barcodes <- input_df %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::mutate(
      ncells = dplyr::n(),
      ncells_to_sample = ifelse(
        ncells >= n_cells_per_class,
        n_cells_per_class,
        ncells
      )
    ) %>%
    dplyr::sample_n(ncells_to_sample, replace = FALSE) %>%
    dplyr::distinct() %>%
    dplyr::pull("barcode")
  input_df <- input_df[selected_barcodes, ]


  # Set parallelization backend
  parallel_trees <- n_trees / n_cpu
  doParallel::registerDoParallel(n_cpu)


  # Fit random forest model
  rf_mod <- foreach::foreach(
    ntree = rep(parallel_trees, n_cpu),
    .combine = randomForest::combine,
    .packages = "randomForest"
  ) %dopar% {
    model <- randomForest::randomForest(
      cluster_id ~ .,
      data = input_df,
      type = "classification",
      importance = TRUE,
      ntree = n_trees,
      verbose = 10
    )
    model
  }


  # Return
  if (return == "oob") {
    oob_acc <- calculate_oob_accuracy(rf_mod)
    oob_acc
  } else if (return == "model"){
    rf_mod
  }
}



find_optimal_resolution <- function(df) {
  df_sorted <- df[order(df$resolution), ]
  score <- c()
  resol <- c()
  oob_accuracies <- df_sorted$oob_accuracy
  resolutions <- df_sorted$resolution
  for (i in 1:(length(resolutions) - 1)) {
    left <- oob_accuracies[1:i]
    right <- oob_accuracies[(i + 1):length(oob_accuracies)]
    var_left <- ifelse(length(left) > 1, var(left), 0)
    var_right <- ifelse(length(right) > 1, var(right), 0)
    diff_means <- mean(left) - mean(right)
    se <- sqrt((var_left / length(left)) + (var_right / length(right)))
    t_score <- diff_means / se
    resol <- c(resol, mean(resolutions[i:(i + 1)]))
    score <- c(score, t_score)
  }
  norm_o <- exp(score) / sum(exp(score))
  df_resol <- data.frame(resol, norm_o)
  colnames(df_resol) <- c("resolution", "optimal_values")
  output <- full_join(df_sorted, df_resol, by = "resolution")
  output
}



plot_split_umap <- function(seurat_obj, var) {
  df <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  df$new_var <- seurat_obj@meta.data[[var]]
  p <- df %>%
    ggplot(aes(UMAP_1, UMAP_2, color = new_var)) +
    geom_point(size = 0.1) +
    facet_wrap(~new_var) +
    theme_classic() +
    theme(legend.position = "none")
  p
}



find_consensus_hvg <- function(seurat_obj, split_by, n_features) {
  seurat_list <- SplitObject(seurat_obj, split.by = split_by)
  seurat_list <- purrr::map(
    seurat_list,
    FindVariableFeatures,
    nfeatures = n_features
  )
  hvg <- purrr::map(seurat_list, VariableFeatures)
  shared_hvg <- Reduce(intersect, hvg)
  shared_hvg
}



plot_annotation_king <- function(seurat_obj, pt_size = 0.25, color_palette) {
  df <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  df$cell_type <- seurat_obj$cell_type
  df$assay <- seurat_obj$assay
  p <- df %>%
    filter(assay == "5P") %>%
    ggplot(aes(UMAP_1, UMAP_2, color = cell_type)) +
    geom_point(size = pt_size) +
    scale_color_manual(values = color_palette[1:length(unique(df$cell_type))]) +
    theme_classic()
  p
}



integrate_object <- function(seurat_obj, n_features, n_dims, group_by_var) {
  seurat_obj <- seurat_obj %>%
    FindVariableFeatures(nfeatures = n_features) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = group_by_var, reduction = "pca", dims = 1:n_dims) %>%
    RunUMAP(reduction = "harmony", dims = 1:n_dims)
  seurat_obj
}



plot_pDNN <- function(seurat_obj, pDNN_var, pt_size = 0.1) {
  df <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  df[[pDNN_var]] <- seurat_obj@meta.data[[pDNN_var]]
  df <- df[!is.na(df[[pDNN_var]]), ]
  p <- ggplot(df, aes_string("UMAP_1", "UMAP_2", color = pDNN_var)) +
    geom_point(size = pt_size) +
    scale_colour_viridis_c(option = "magma") +
    theme_classic()
  p
}



seurat2shiny = function(
  object                         ,
  assay     = object@active.assay,
  slot      = "data"             ,
  reduction = "umap"
) {

  # Input check.
  if ( ! is(object, "Seurat") )
    stop("'object' is not a Seurat object.");
  
  if ( ! assay %in% Seurat::Assays(object) )
    stop("'assay' not in the Seurat object's available assays.");
  
  if ( ! slot %in% c("counts", "data", "scale.data") )
    stop("'slot' not in the Seurat object's available slots.");

  
  # Extract expression data.
  expression = Seurat::GetAssayData(object = object, slot = slot, assay = assay);
  
  
  # Extract 2D coordinates.
  embeds <- as.data.frame(object@reductions[[reduction]]@cell.embeddings);
  names(embeds) <- c("coord_x", "coord_y");

  
  # Join metadata with coordinates.
  metadata <- object@meta.data;
  if (all(rownames(metadata) == rownames(embeds))) {
    metadata$coord_x <- embeds$coord_x
    metadata$coord_y <- embeds$coord_y
    metadata$barcode <- rownames(metadata)
  }

  
  if ( ! identical( as.character(metadata$barcode), colnames(expression) ) )
    warning("Cells in metadata and expression matrix do not match.");
  

  invisible(
    list(metadata = metadata, expression = expression)
  );
}


find_proportions_df <- function(seurat_obj, x, fill) {
  df <- seurat_obj@meta.data %>%
    select(x, fill) %>%
    group_by(.data[[x]], .data[[fill]]) %>%
    summarise(n_cells = n()) %>%
    ungroup() %>%
    group_by(.data[[x]]) %>%
    mutate(n_cells_total = sum(n_cells)) %>%
    ungroup() %>%
    mutate(percentage_cells = round(n_cells / n_cells_total * 100, 3))
  df
}


plot_stacked_barplot <- function(df, x, fill, colors) {
  p <- df %>%
    ggplot(aes_string(x, "percentage_cells", fill = fill)) +
    geom_col() +
    ggtitle("Age group") +
    labs(x = "", y = "Percentage of Cells (%)", fill = "") + 
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5, face = "bold")
    )
  p
}


find_markers_subclusters <- function(seurat_obj, var, pattern) {
  clusters_interest <- seurat_obj@meta.data[, var] %>%
    unique() %>%
    str_subset(pattern) %>%
    sort()
  markers_interest <- purrr::map(clusters_interest, function(x) {
    group_1 <- clusters_interest[which(clusters_interest == x)]
    group_2 <- clusters_interest[which(clusters_interest != x)]
    df <- FindMarkers(
      seurat_obj,
      ident.1 = group_1,
      ident.2 = group_2,
      only.pos = TRUE,
      logfc.threshold = 0.5,
      verbose = TRUE
    )
    df <- df %>%
      rownames_to_column(var = "gene") %>%
      arrange(desc(avg_log2FC))
    df
  })
  names(markers_interest) <- clusters_interest
  markers_interest
}


updateAnnotation <- function(seurat_obj,
                            refAnnotation = "20220215",
                            newAnnotation = "20220619") {
  dict_20220619 <- c(
    "GC-Tfh-0X40" = "GC-Tfh-OX40",
    "non-GC-Tf-regs" = "Eff-Tregs-IL32",
    "GC-Tf-regs" = "Tfr",
    "IFN CD8 T" = "IFN+ CD8 T",
    "CXCR6+ RM CD8 T" = "RM CD8 activated T",
    "IL7R MMP12 macrophages" = "MMP Slancytes",
    "C1Q HLA macrophages" = "C1Q Slancytes",
    "SELENOP FUCA1 PTGDS macrophages" = "SELENOP Slancytes",
    "ITGAX ZEB2 macrophages" = "ITGAX Slancytes",
    "Mast cells" = "Mast",
    "LZ FDC" = "FDC",
    "DZ FDC" = "COL27A1+ FDC",
    "DZ_Sphase" = "DZ early Sphase",
    "DZ_Sphase_HistoneHigh" = "DZ late Sphase",
    "DZ_G2M_HistoneHigh" = "DZ early G2Mphase",
    "DZ_G2M_CCNBHigh"= "DZ late G2Mphase",
    "DZ-cell cycle exit" = "DZ cell cycle exit",
    "DZ-nonproliferative" = "DZ cell cycle exit",
    "DZ-nonproliferative_FOXP1hi"= "DZ non proliferative",
    "DZ/LZ" = "DZ_LZ transition",
    "DZ/LZ" = "DZ_LZ transition",
    "LZ" = "LZ",
    "LZ-BCL2A1 neg"= "LZ",
    "LZ-DZ-re-entry early commitment" = "LZ_DZ reentry commitment",
    "LZ-proliferative_BCL2A1pos" = "LZ proliferative",
    "LZ-proliferative_BCL2A1neg" = "LZ_DZ transition",
    "MBC-like_nonproli" = "Precursor MBCs",
    "MBC-like_proli1" = "Precursor MBCs",
    "MBC-like_proli2"= "Reactivated proliferative MBCs",
    "MBC-like_proli3" = "Reactivated proliferative MBCs",
    "MBC-like_FCRL4+"= "Reactivated proliferative MBCs",
    "PC-precursors" = "PC committed Light Zone GCBC",
    "class switch MBC" = "csMBC",
    "Neutrophil Granulocytes" = "Neutrophils",
    
    
    # Update post-revision
    "MAIT" = "MAIT/Vδ2+ γδ T",
    "TCRVδ+ gd T" = "non-Vδ2+ γδ T",
    "CD56+ gd T" = "ZNF683+ CD8 T",
    "Nksig CD8 T" = "EM CD8 T"
  )
  oldAnnot <- paste("annotation", refAnnotation, sep = "_")
  newAnnot <- paste("annotation", newAnnotation, sep = "_")
  seurat_obj@meta.data[[newAnnot]] <- as.character(seurat_obj@meta.data[[oldAnnot]])
  dictSub <- dict_20220619[names(dict_20220619) %in% unique(seurat_obj@meta.data[[oldAnnot]])]
  for (cellType in names(dictSub)) {
    seurat_obj@meta.data[[newAnnot]][seurat_obj@meta.data[[newAnnot]] == cellType] <- dictSub[cellType]
  }
  seurat_obj
}