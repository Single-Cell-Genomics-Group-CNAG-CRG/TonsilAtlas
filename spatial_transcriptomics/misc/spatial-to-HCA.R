## Load libraries
library(Seurat)
library(dplyr)
library(tibble)
library(glue)
library(here)

## Set parameters
set.seed(123)
source(here("misc/paths.R"))

"misc/{plt_dir}" %>%
    glue() %>%
    here() %>%
    dir.create(path = .,
        showWarnings = FALSE,
        recursive = TRUE)

"misc/{robj_dir}" %>%
    glue() %>%
    here() %>%
    dir.create(path = .,
        showWarnings = FALSE,
        recursive = TRUE)

## Load data
# The data used in this Rmarkdown document comes from **03-clustering_integration.Rmd** where the data was integrated.
merged_se <- "{anot}/{robj_dir}/integrated_spatial_annot.rds" %>%
    glue() %>%
    here() %>%
    readRDS(file = .)

meta <- "01-spaceranger/data/tonsil_atlas_donor_metadata.csv" %>%
    glue() %>%
    here() %>%
    readr::read_csv()

sample <- "01-spaceranger/data/sample_id.txt" %>%
    glue() %>%
    here() %>%
    readr::read_csv()
    

# Extract UMAP coordinates
umap_coord <- Seurat::Reductions(merged_se, "umap")@cell.embeddings %>%
    data.frame %>%
    rename(
        UMAP_1_20220215 = UMAP_1,
        UMAP_2_20220215 = UMAP_2) %>%
    tibble::rownames_to_column("barcode")

# Cell cycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_se <- CellCycleScoring(
    merged_se,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE)

# Columns to kep
keep <- c("barcode", "donor_id", "gem_id", "library_name", "assay", "sex", "age",
    "age_group", "hospital", "nCount_Spatial", "nFeature_Spatial", "pct_mt",
    "pct_ribosomal","S.Score", "G2M.Score", "Phase",
    "annotation_20220215", "UMAP_1_20220215", "UMAP_2_20220215", "index",
    "slide", "area")

merged_se@meta.data <- merged_se@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    select(-c("hospital", "sex", "age")) %>%
    left_join(meta, by = "donor_id") %>%
    left_join(sample, by = c("donor_id", "gem_id")) %>%
    left_join(umap_coord, by = "barcode") %>% 
    rename(
        annotation_20220215 = annotation.general,
        pct_mt = percent.mito,
        pct_ribosomal = percent.ribo,
        library_name = library_id,
        assay = type) %>%
    select(all_of(keep)) %>%
    tibble::column_to_rownames("barcode")

"misc/{robj_dir}/20220215_tonsil_atlas_spatial_seurat_obj.rds.rds" %>%
    glue() %>%
    here() %>%
    saveRDS(object = merged_se, file = .)
