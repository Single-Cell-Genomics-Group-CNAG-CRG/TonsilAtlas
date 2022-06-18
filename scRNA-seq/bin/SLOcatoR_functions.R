# This script contains functions to project query single-cell transcriptomes
# against the tonsil atlas referece. In particular, it allows to:
# 1. Integrate query and reference datasets "de novo" using Harmony.
# 2. Transfers the labels of the reference to the query using KNN classifier.
# 3. Transfers the UMAP coordinates via KNN regression to conserve the embeddings


find_assay_specific_features <- function(seurat_obj,
                                         assay_var = "assay",
                                         n_features = 5000) {
  seurat_list <- SplitObject(seurat_obj, split.by = assay_var)
  seurat_list <- purrr::map(
    seurat_list,
    FindVariableFeatures,
    nfeatures = n_features
  )
  hvg <- purrr::map(seurat_list, VariableFeatures)
  shared_hvg <- Reduce(intersect, hvg)
  shared_hvg
}


integrate_assays <- function(seurat_obj,
                             assay_specific = TRUE,
                             assay_var = "assay",
                             shared_hvg,
                             n_dim = 30
) {
  if (assay_specific) {
    seurat_obj <- seurat_obj %>%
      ScaleData(features = shared_hvg) %>%
      RunPCA(features = shared_hvg) %>%
      RunHarmony(group.by.vars = assay_var, reduction = "pca", dims = 1:n_dim)
  } else {
    seurat_obj <- seurat_obj %>%
      ScaleData() %>%
      RunPCA()
  }
  
  seurat_obj
}


split_training_and_test_sets <- function(seurat_obj,
                                         split_var,
                                         referece_label,
                                         query_label,
                                         reduction = "harmony",
                                         n_dims = 30) {
  # Training and test sets are composed by the reference and query cells
  cells <- colnames(seurat_obj)
  reference_cells <- cells[seurat_obj@meta.data[[split_var]] == referece_label]
  query_cells <- cells[seurat_obj@meta.data[[split_var]] == query_label]
  
  
  # Batch-corrected principal components are used as explanatory variables
  cell_embeddings <- Embeddings(seurat_obj, reduction = "harmony")
  training_set <- cell_embeddings[reference_cells, 1:n_dims]
  test_set <- cell_embeddings[query_cells, 1:n_dims]
  list(training_set = training_set, test_set = test_set)
}


find_optimal_k <- function(seurat_obj,
                           training_set,
                           response_var,
                           ks = c(2, 4, 6, 8, 16, 32, 64, 128, 256),
                           return_plot = TRUE) {
  indices <- sample(
    1:nrow(training_set),
    size = nrow(training_set) * 0.7,
    replace = FALSE
  ) 
  train_loan  <- training_set[indices, ] # 70% training data
  test_loan <- training_set[-indices, ] # remaining 30% test data
  train_labels <- seurat_obj@meta.data[rownames(train_loan), ][[response_var]]
  test_labels <- seurat_obj@meta.data[rownames(test_loan), ][[response_var]]
  
  k_optm <- c()
  k_values <- c()
  
  for (i in ks) {
    print(i)
    knn_mod <- knn(
      train = train_loan,
      test = test_loan,
      cl = train_labels,
      k = i
    )
    k_optm <- c(k_optm, mean(test_labels == knn_mod) * 100)
    k_values <- c(k_values, i)
  }
  
  k_df <- data.frame(k = k_values, accuracy = k_optm)
  
  if (return_plot) {
    p <- ggplot(k_df, aes(k, accuracy)) +
      geom_line() +
      geom_point() +
      theme_bw()
    return(list(df = k_df, plot = p))
    
  } else {
    return(k_df)
  }
}


transfer_label <- function(seurat_obj,
                           training_set,
                           test_set,
                           response_var,
                           k) {
  training_labels <- seurat_obj@meta.data[rownames(training_set), response_var]
  knn_mod <- knn(
    train = training_set,
    test = test_set,
    cl = training_labels,
    k = k,
    prob = TRUE
  )
  test_labels_df <- data.frame(
    query_cells = rownames(test_set),
    annotation = as.character(knn_mod),
    annotation_prob = attr(knn_mod, "prob")
  )
  test_labels_df
}


transfer_umap_coords <- function(seurat_obj,
                                 training_set,
                                 test_set,
                                 umap1_var,
                                 umap2_var,
                                 k) {
  train_umap_1 <- seurat_obj@meta.data[rownames(training_set), umap1_var]
  train_umap_2 <- seurat_obj@meta.data[rownames(training_set), umap2_var]
  knn_mod_umap_1 <- knnreg(x = training_set, y = train_umap_1, k = k)
  knn_mod_umap_2 <- knnreg(x = training_set, y = train_umap_2, k = k)
  umap_1_pred <- predict(knn_mod_umap_1, newdata = as.data.frame(test_set))
  umap_2_pred <- predict(knn_mod_umap_2, newdata = as.data.frame(test_set))
  umap_pred_df <- data.frame(
    query_cells = rownames(test_set),
    UMAP1 = umap_1_pred,
    UMAP2 = umap_2_pred
  )
  umap_pred_df
}


integrate_assays_atac <- function(seurat_obj,
                                  n_dim = 40,
                                  assay_use,
                                  group_by_vars,
                                  reduction) {
  seurat_obj <- seurat_obj %>%
    RunTFIDF() %>% 
    FindTopFeatures(min.cutoff = "q0") %>%
    RunSVD() %>%
    RunHarmony(
      group.by.vars = group_by_vars, 
      reduction = reduction,
      dims = 2:n_dim,
      assay.use = assay_use,
      project.dim = FALSE
    )
  seurat_obj
}


install_tonsil <- function(path_to_save) {
  dataset_url <- "https://zenodo.org/record/5186413/files/TICAtlas_metadata.csv"
  command <- glue::glue("wget -O {path_to_save} {dataset_url}")
  system(command)
}
