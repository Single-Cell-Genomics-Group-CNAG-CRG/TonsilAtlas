library(Seurat)

source(here::here("misc/paths.R"))
source(here::here("utils/bin.R"))

# Load data
merged_se <- "{clust}/{robj_dir}/integrated_spatial.rds" %>%
  glue::glue() %>%
  here::here() %>%
  readRDS(file = .)

# Add metadata
merged_se[["annotation"]] <- dplyr::case_when(
  merged_se@meta.data$Spatial_snn_res.0.3 == 0 ~ "Inter-follicular zone 1",
  merged_se@meta.data$Spatial_snn_res.0.3 == 1 ~ "T cell zone",
  merged_se@meta.data$Spatial_snn_res.0.3 == 2 ~ "Follicle",
  merged_se@meta.data$Spatial_snn_res.0.3 == 3 ~ "Epithelial 1",
  merged_se@meta.data$Spatial_snn_res.0.3 == 4 ~ "Follicle Proliferating",
  merged_se@meta.data$Spatial_snn_res.0.3 == 5 ~ "Epithelial 2",
  merged_se@meta.data$Spatial_snn_res.0.3 == 6 ~ "Inter-follicular zone 2",
  merged_se@meta.data$Spatial_snn_res.0.3 == 7 ~ "Muscle",
)


# Save individual objects
lapply(Seurat::Images(merged_se), function(img) {
  
  # Subset Seurat object
  se_sub <- subset(merged_se, subset = gem_id == "esvq52_nluss5")
  se_sub
  se_sub@images <- se_sub@images[Seurat::Images(se_sub) == "esvq52_nluss5"]
  
  # extract spatial coordinates
  embeds <- data.frame(se_sub@images[[img]]@coordinates[, c("imagerow", "imagecol")])
  colnames(embeds) <- c("coord_y", "coord_x");
  # Inverse coord_y
  embeds$coord_y <- - embeds$coord_y
  
  # Join metadata with coordinates.
  metadata <- se_sub@meta.data;
  
  metadata <- merge(x = metadata, y = embeds, by = "row.names");
  names(metadata)[1] <-  "barcode"; # names(metadata)[names(metadata) == "Row.names"] = "barcode";
  rownames(metadata) <- metadata$barcode
  metadata$barcode <- as.character(metadata$barcode)
  # Reset barcode order which is shuffled in merge
  metadata <- metadata[rownames(se_sub@meta.data), ]
  
  # Save Metadata
  "misc/metadata_{img}.rds" %>%
    glue::glue() %>%
    here::here() %>%
    saveRDS(object = metadata, file = .)
  
  
  })
