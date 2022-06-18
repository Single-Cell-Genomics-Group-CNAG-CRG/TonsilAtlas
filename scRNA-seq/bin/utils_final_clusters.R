# This script defines the final paths for each object needed to plot the UMAP
# in figure 1

# Paths
path_to_final_myeloid <- here("scRNA-seq/results/R_objects/level_5/myeloid/myeloid_annotated_level_5.rds")
path_to_multiome_myeloid <- here("scRNA-seq/results/R_objects/level_2/myeloid/myeloid_clustered_filtered_level_2.rds")

path_to_final_epithelial <- here("scRNA-seq/results/R_objects/level_4/epithelial/epithelial_annotated_level_4.rds")
path_to_multiome_epithelial <- here("scRNA-seq/results/R_objects/level_2/epithelial/epithelial_clustered_filtered_level_2.rds")
path_to_fdc_with_epithelial <- here("scRNA-seq/results/R_objects/level_4/FDC/epithelial_level_4.rds")

path_to_final_pdc <- here("scRNA-seq/results/R_objects/level_4/PDC/PDC_annotated_level_4.rds")
path_to_multiome_pdc <- here("scRNA-seq/results/R_objects/level_2/PDC/PDC_clustered_filtered_level_2.rds")

path_to_final_fdc <- here("scRNA-seq/results/R_objects/level_5/FDC/FDC_all_level_5.rds")
path_to_multiome_fdc <- here("scRNA-seq/results/R_objects/level_2/FDC/FDC_clustered_filtered_level_2.rds")
path_to_fdc_proliferative <- here("scRNA-seq/results/R_objects/level_4/FDC/FDC_proliferative_level_4.rds")

path_to_preB <- here("scRNA-seq/results/R_objects/level_2/preBC/preBC_subsetted_level_2.rds")
path_to_preT <- here("scRNA-seq/results/R_objects/level_2/preTC/preTC_subsetted_level_2.rds")

path_to_cd8_T <- here("scRNA-seq/results/R_objects/level_5/Cytotoxic/paper/CD8_T_level_5_annotated_level_5.rds")
path_to_ilc_or_nk <- here("scRNA-seq/results/R_objects/level_5/Cytotoxic/ILC_NK/ILC_NK_annotated_level_5.rds")

path_to_cd4 <- here("scRNA-seq/results/R_objects/level_5/CD4_T/paper/CD4_T_subseted_annotated_level_5.rds")
path_to_proliferative_cd4_t <- here("scRNA-seq/results/R_objects/level_5/CD4_T/CD4_T_proliferative_subsetted_level_5.rds")

path_to_nbc_mbc <- here("scRNA-seq/results/R_objects/final_clusters/NBC_MBC_seu_obj_level_5_delta.rds")
path_to_gcbc <- here("scRNA-seq/results/R_objects/final_clusters/GCBC_seu_obj_level_5_beta.rds")
path_to_pc <- here("scRNA-seq/results/R_objects/final_clusters/PC_seu_obj_level_5_eta.rds")
path_to_pc_lvl_3 <- here("scRNA-seq/results/R_objects/level_3/PC/PC_clustered_level_3.rds")
path_to_save <- here("scRNA-seq/results/R_objects/final_clusters/dataframe_with_all_cells_20220215.rds")

path_to_th <- here("scRNA-seq/results/R_objects/level_5/CD4_T/paper/CD4_T_Th_level_5_annotated.rds")

# Paths to save harmonized objects
path_to_save_myeloid <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_myeloid_seurat_obj.rds")
path_to_save_epithelial <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_epithelial_seurat_obj.rds")
path_to_save_pdc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_PDC_seurat_obj.rds")
path_to_save_fdc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_FDC_seurat_obj.rds")
path_to_save_preB <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_preB_seurat_obj.rds")
path_to_save_preT <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_preT_seurat_obj.rds")
path_to_save_cd8 <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_CD8_T_seurat_obj.rds")
path_to_save_ilc_nk <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_ILC_NK_seurat_obj.rds")
path_to_save_cd4 <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_CD4_T_seurat_obj.rds")
path_to_save_nbc_mbc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_NBC_MBC_seurat_obj.rds")
path_to_save_gcbc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_GCBC_seurat_obj.rds")
path_to_save_pc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_PC_seurat_obj.rds")
path_to_save_th <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_Th_seurat_obj.rds")

path_to_save_tonsil <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_tonsil_atlas_rna_seurat_obj.rds")
path_to_save_tonsil_cite <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_tonsil_atlas_cite_seurat_obj.rds")
path_to_save_tonsil_atac <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_tonsil_atlas_atac_seurat_obj.rds")

path_to_spatial <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_tonsil_atlas_spatial_seurat_obj.rds")

path_to_save_umap_df_csv <- here("data/raw_data_figures/all_cells_rna_multiome_to_create_supp_table_1.csv")
path_to_save_umap_df_rds <- here("scRNA-seq/3-clustering/final_clusters/tmp/umap_df_fig1_rna_multiome.rds")
path_to_save_df_multi_myeloid <- here("scRNA-seq/3-clustering/final_clusters/tmp/umap_df_fig1_rna_multiome_myeloid_annotated.rds")
path_to_save_cite_seq_df <- here("scRNA-seq/3-clustering/final_clusters/tmp/umap_df_fig1_rna_multiome_myeloid_annotated_and_cite_seqs.rds")
path_to_save_atac_seq_df <- here::here("scRNA-seq/3-clustering/final_clusters/tmp/umap_df_atac_multiome.rds")


path_to_tonsil_level_1 <- here("scRNA-seq/results/R_objects/tonsil_rna_integrated_annotated_level_1.rds")
path_to_cite_seq <- here("CITE-seq/results/seurat_object_cite_seq_seurat_wnn.rds")
path_to_atac_seq <- here("scATAC-seq/results/R_objects/8.3.tonsil_peakcalling_annotation_level1_signature.rds")


path_to_df_umap_fig1 <- here("data/raw_data_figures/figure_1_all_technologies_umaps.csv")
path_to_qc_metrics_rna_multi_cite <- here("data/raw_data_figures/supp_table_qc_metrics_rna_multiome_cite.csv")
path_to_qc_metrics_atac_multi <- here("data/raw_data_figures/supp_table_qc_metrics_atac_multiome.csv")


path_to_save_umaps_fig1 <- here("results/paper/figures/figure_1_main_umaps.pdf")

path_to_save_qc_table <- here("results/paper/tables/supplementary_table_qc_metrics.xlsx")


path_to_cite_cd4_not_harmonized <- here("CITE-seq/results/R_objects/final_clusters/20220215/CD4_tonsil_cite_seq_annotated_tcr_analysis_updated_obj.rds")
path_to_cite_cd8_not_harmonized <- here("scRNA-seq/results/R_objects/level_5/Cytotoxic/paper/CD8_T_level_5_integrated_all_CITEseq.rds")
path_to_save_cite_cd4 <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_CD4_T_seurat_obj_cite.rds")
path_to_save_cite_cytotoxic <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_cytotoxic_seurat_obj_cite.rds")
path_to_atac_cd4_not_harmonized <- here("data/raw_data_figures/scATAC_CD4_T/files_plots/CD4T_integration_peak_calling_level_5.rds")
path_to_save_atac_cd4 <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_CD4_T_seurat_obj_atac.rds")
# path_to_atac_cytotoxic_not_harmonized <- here("scATAC-seq/results/R_objects/05.Cytotoxic_chromVar_CISBP_level_4.rds")
path_to_atac_cytotoxic_not_harmonized <- here("scATAC-seq/results/R_objects/05.Cytotoxic_chromVar_JASPAR_level_4.rds")
path_to_save_atac_cytotoxic <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_cytotoxic_seurat_obj_atac.rds")
path_to_atac_nbc_mbc_not_harmonized <- here("scATAC-seq/results/R_objects/Targetted_analysis/NBC_MBC/NBC_MBC_integrated_level_3.rds")
path_to_atac_gcbc_not_harmonized <- here("scATAC-seq/results/R_objects/Targetted_analysis/GCBC/GCBC_chromVar_CISBP_level_4.rds")
path_to_atac_pc_not_harmonized <- here("scATAC-seq/results/R_objects/Targetted_analysis/PC/05.PC_chromVar_CISBP_level_5.rds")
path_to_save_atac_nbc_mbc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_NBC_MBC_seurat_obj_atac.rds")
path_to_save_atac_gcbc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_GCBC_seurat_obj_atac.rds")
path_to_save_atac_pc <- here("scRNA-seq/results/R_objects/final_clusters/20220215/20220215_PC_seurat_obj_atac.rds")


# Variables
selected_cols <- c("barcode", "donor_id", "gem_id", "library_name","assay",
                   "sex", "age", "age_group", "hospital", "UMAP_1_level_1", "UMAP_2_level_1",
                   "nCount_RNA", "nFeature_RNA", "pct_mt", "pct_ribosomal", "is_hashed",
                   "pDNN_hashing", "pDNN_scrublet", "pDNN_union",
                   "scrublet_doublet_scores", "scrublet_predicted_doublet",
                   "S.Score", "G2M.Score", "Phase", "CC.Difference",
                   "annotation_level_1","annotation_figure_1",
                   "annotation_20220215", "UMAP_1_20220215", "UMAP_2_20220215")
selected_cols_cite <- c("barcode", "donor_id", "subproject", "gem_id","assay", "sex",
                        "age", "age_group", "hospital", "nCount_RNA", "nCount_ADT",
                        "nFeature_ADT", "nFeature_RNA", "pct_mt", "pct_ribosomal",
                        "scrublet_doublet_scores", "scrublet_predicted_doublet",
                        "S.Score", "G2M.Score", "Phase", "UMAP_1_20220215",
                        "UMAP_2_20220215", "bcr_flag", "tcr_flag", "RNA.weight",                    
                        "ADT.weight")
selected_cols_atac <- c("barcode", "donor_id", "gem_id", "library_name", "assay", "sex", "age", "age_group",
                        "hospital", "is_facs_sorted", "UMAP_1_level_1", "UMAP_2_level_1",
                        "nCount_ATAC", "nFeature_ATAC", "nCount_peaks_macs", "nFeature_peaks_macs",
                        "TSS_fragments", "TSS.enrichment", "TSS.percentile", "high.tss", "DNase_sensitive_region_fragments",
                        "enhancer_region_fragments", "promoter_region_fragments", "on_target_fragments",
                        "blacklist_region_fragments", "blacklist_ratio", "peak_region_fragments", "peak_region_cutsites",
                        "pct_reads_in_peaks", "nucleosome_signal", "nucleosome_percentile", "duplicate",
                        "scrublet_doublet_scores_atac", "scrublet_predicted_doublet_atac", "annotation_level_1",
                        "annotation_figure_1",  "annotation_prob","UMAP_1_20220215", "UMAP_2_20220215")




color_palette <-  c("#1CFFCE", "#90AD1C", "#C075A6", "#85660D", "#5A5156", "#AA0DFE",   
                    "#F8A19F", "#F7E1A0", "#1C8356", "#FEAF16", "#822E1C", "#C4451C",   
                    "#1CBE4F", "#325A9B", "#F6222E", "#FE00FA", "#FBE426", "#16FF32", 
                    "black",   "#3283FE", "#B00068", "#DEA0FD", "#B10DA1", "#E4E1E3",   
                    "#90AD1C", "#FE00FA", "#85660D", "#3B00FB", "#822E1C", "coral2", 
                    "#1CFFCE", "#1CBE4F", "#3283FE", "#FBE426", "#F7E1A0", "#325A9B",   
                    "#2ED9FF", "#B5EFB5", "#5A5156", "#DEA0FD", "#FEAF16", "#683B79",   
                    "#B10DA1", "#1C7F93", "#F8A19F", "dark orange", "#FEAF16", "#FBE426",  
                    "Brown")


cols_fig1 <- c(
  "NBC" = "#dcf0f4",
  "Activated NBC" = "#7bc6d6", 
  "GCBC" = "#398D9F",
  "MBC" = "#025566", 
  "PC" = "#032872", 
  
  "CD4 T" = "#a50000",
  "Naive CD4 T" = "gray50", 
  "Naive CD8 T" = "#dfa6b8", 
  "CD8 T" = "#bf4d72", 
  "NK" = "#67253a", 
  "DN" = "#a55a88", 
  "ILC" = "#8c324f",
  "cycling T"  = "#cc0000",
  
  "preB/T" = "#4363d8", 
  
  
  "DC" = "#b6dec3", 
  "Mono/Macro" = "#49ae6b", 
  "Mast" = "#27723c", 
  "Granulocytes" = "#00cd00", 
  
  "FDC" = "#d5ff00", 
  "cycling FDC" = "#aacc00", 
  
  "epithelial" = "#ffdc9f", 
  
  "PDC" = "#f032e6", 
  
  "cycling myeloid" = "black"
)

# Old paths (when run locally)
# path_to_final_myeloid <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/myeloid/myeloid_annotated_level_5.rds"
# path_to_multiome_myeloid <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_2/myeloid/myeloid_clustered_filtered_level_2.rds"
# 
# path_to_final_epithelial <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/epithelial/epithelial_annotated_level_4.rds"
# path_to_multiome_epithelial <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_2/epithelial/epithelial_clustered_filtered_level_2.rds"
# path_to_fdc_with_epithelial <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/FDC/epithelial_level_4.rds"
# 
# path_to_final_pdc <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/PDC/PDC_annotated_level_4.rds"
# path_to_multiome_pdc <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_2/PDC/PDC_clustered_filtered_level_2.rds"
# 
# path_to_final_fdc <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/FDC/FDC_all_level_4.rds"
# path_to_multiome_fdc <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_2/FDC/FDC_clustered_filtered_level_2.rds"
# path_to_fdc_proliferative <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/FDC/FDC_proliferative_level_4.rds"
# 
# path_to_preB <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_2/preBC/preBC_subsetted_level_2.rds"
# 
# path_to_preT <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_2/preTC/preTC_subsetted_level_2.rds"
# 
# path_to_cd8_T <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/Cytotoxic/CD8_T/CD8_T_clustered_level_5.rds"
# path_to_ilc_or_nk <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/Cytotoxic/ILC_NK/ILC_NK_clustered_level_5.rds"
# path_to_tim3_pos_T <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/Cytotoxic/TIM3_DN/TIM3_DN_clustered_level_4.rds"
# 
# path_to_cd4 <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/CD4_T/CD4_T_subseted_integrated_level_5_2.rds"
# path_to_proliferative_cd4_t <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/CD4_T/CD4_T_proliferative_subsetted_level_5.rds"
# path_to_naive_cd8 <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/CD4_T/naive_CD8_T_subsetted_level_4.rds"
# 
# path_to_b_cells <- "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_4/B_cells/Bcells_level4.rds"



