# A periodic table of human tonsillar cells

Palatine tonsils are secondary lymphoid organs representing the first line of immunological defense against inhaled or ingested pathogens. Here, we provide a comprehensive census of cell types forming the human tonsil by applying single cell transcriptome, epigenome, proteome and adaptive immune repertoire as well as spatial transcriptomics, resulting in an atlas of >357,000 cells. We provide a glossary of 137 annotated cell types and states, and disentangle gene regulatory mechanisms that drive cells through specialized lineage trajectories. Exemplarily, we stratify multiple tonsil-resident myeloid slancyte subtypes, establish a distant BCL6 superenhancer as locally active in both follicle-associated B and T cells, and describe SIX5 as a potentially novel transcriptional regulator of plasma cell maturation. Further, our atlas acts as a reference map to pinpoint alterations observed in disease, shown here discovering immune-phenotype plasticity in tumoral cells and microenvironment shifts of mantle cell lymphomas (MCL). To facilitate such reference-based analysis, we develop HCATonsilData and SLOcatoR, a computational framework that provides programmatic and modular access to our dataset; and allows the straightforward annotation of future single-cell profiles from secondary lymphoid organs.


This repository contains all the scripts, notebooks and reports to reproduce all analysis from our recent preprint titled "A periodic table of tonsillar cells" (TODO: update title). Here, we describe how to access the data, document the most important packages and versions used, and explain how to navigate the directories and files in this repository.

![](data/TonsilAtlasPic.png.png)


## Data

The data has been deposited in five levels of organization, from raw to processed data:


* Level 1: raw data. All fastq files for all data modalities have been deposited at the European Genome Archive (EGA) under accession id XXX.
* Level 2: matrices. All data modalities correspond to different technologies from 10X Genomics. As such, they were mapped with different flavors of CellRanger (CR). The most important files in the “outs” folder of every CR run (including all matrices) have been deposited in Zenodo (TODO: include link).
* Level 3: Seurat Objects. All data was analyzed within the Seurat ecosystem. We have archived in [Zenodo](https://doi.org/10.5281/zenodo.6340174) all Seurat Objects that contain the raw and processed counts, dimensionality reductions (PCA, Harmony, UMAP), and metadata needed to reproduce all figures from this manuscript.
* Level 4: to allow for programmatic and modular access to the whole tonsil atlas dataset, we developed HCATonsilData, available on [GitHub](https://github.com/massonix/HCATonsilData). HCATonsilData provides a vignette which documents how to navigate and understand the data. In addition, we will periodically update the annotations as we refine it with suggestions from the  community.
* Level 5: interactive mode. We provide different iSEE instances to browse the data interactively: [http://shiny.imbei.uni-mainz.de:3838/iSEE_TonsilDataAtlas/](http://shiny.imbei.uni-mainz.de:3838/iSEE_TonsilDataAtlas/).


We refer to the READMEs in these repositories for a full explanation of the dataset.


## Package versions

* [Seurat 3.2.0 and 4.1.0](https://satijalab.org/seurat/)
* [Signac 1.1.0](https://github.com/timoast/signac)
* [harmony 1.0](https://github.com/immunogenomics/harmony)
* [lisi 1.0](https://github.com/immunogenomics/LISI)
* [scrublet 0.2.1](https://github.com/swolock/scrublet)
* [UCell 1.99.7](https://github.com/carmonalab/UCell)
* [clusterProfiler 4.3.4](https://yulab-smu.top/biomedical-knowledge-mining-book/)
* [pySCENIC 0.10.3](https://scenic.aertslab.org/)
* [CellPhoneDB 3](https://www.cellphonedb.org/)
* [chromVar 1.1.0](https://bioconductor.org/packages/release/bioc/html/chromVAR.html)
* [JASPAR2020 0.99.10](http://bioconductor.org/packages/release/data/annotation/html/JASPAR2020.html)
* [chromVARmotifs 0.2.0](https://github.com/GreenleafLab/chromVARmotifs)
* [Vireo 0.5.0](https://github.com/single-cell-genetics/vireo)
* [Scirpy 0.7.0](https://scverse.org/scirpy/latest/tutorials.html)
* [SPOTlight 0.1.7](https://bioconductor.org/packages/release/bioc/html/SPOTlight.html)
* [Rmagic 2.0.3](https://github.com/KrishnaswamyLab/MAGIC)
* [SPATA2 0.1.0](https://github.com/theMILOlab/SPATA2)

You can check the versions of other packages at the "Session Information" section of each html report. To visualize one of the html reports online, you can copy&paste the URL of the report directly into the [HTML GitHub viewer](https://htmlpreview.github.io/).


# File system and name scheme

Although each technology requires specific analysis, they also share a similar pre-processing pipeline.
We have strived to harmonize these pipelines into similar naming schemes so that it is easy for users
to navigate this repo. Likewise, we have tried to code in a shared style. These are the most important steps:

* 1-cellranger_mapping: scripts used to run cellranger in our cluster. It also contains QC metrics for different sequencing runs.
* 2-QC: quaity control for the sequencing and mapping of raw data, filtering of poor-quality cells and genes, normalization, doublet detection, and batch effect correction.
* 3-clustering: we followed a top-down clustering approach (see methods of our manuscript). Thus, the clustering is organized by levels, in which we move from general cellular compartments to granular cell types and states in a hierarchical and recursive fashion. In this notebooks we have also included the annotation, which was established in collaboration with the annotation team.

Other more focused analysis include:

* [Gene regulatory inference using pySCENIC](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/tree/main/scRNA-seq/gene_regulatory_networks/pyScenic)
* [Gene Set Enrichment analysis for different slancyte subsets](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/blob/main/scRNA-seq/gsea/gsea_myeloid.Rmd)
* [Cell-cell interactions with cellPhoneDB](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/tree/main/scRNA-seq/cell_to_cell_networks)
* [Label transfer from RNA to ATAC](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/blob/main/scATAC-seq/3-Label_transfering/01-scATAC_annotation_KNN.Rmd)
* [Motif analysis with chromVAR](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/tree/main/scATAC-seq/5-motif_analysis)
* [Clonotype analysis with scirpy](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/blob/main/CITE-seq/08-sub_clustering/CD4_T/05-scirpy_tcr.ipynb)
* [Spot deconvolution with SPOTlight](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/tree/main/spatial_transcriptomics/05-sc_map)


In addition, the ["figures_and_scripts"](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/tree/main/scRNA-seq/figures_and_tables) folder contains the scripts used to generate most of the figures in the manuscript. Finally, the ["bin"](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/tree/main/scRNA-seq/bin) folder contains functions and utilities used throughout many scripts.

# Getting the code

You can download a copy of all the files in this repository by cloning the git repository:

```{bash}
git clone https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas.git
```

# Relevant literature

* [Single-cell analysis of human B cell maturation predicts how antibody class switching shapes selection dynamics](https://www.science.org/doi/10.1126/sciimmunol.abe6291?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
* [Dynamics of B cells in germinal centres](https://www.nature.com/articles/nri3804)
* [The generation of antibody-secreting plasma cells](https://www.nature.com/articles/nri3795)
* [Integrated single-cell transcriptomics and epigenomics reveals strong germinal center–associated etiology of autoimmune risk loci](https://www.science.org/doi/10.1126/sciimmunol.abh3768?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
* [Characterization of human FDCs reveals regulation of T cells and antigen presentation to B cells](https://rupress.org/jem/article/218/10/e20210790/212590/Characterization-of-human-FDCs-reveals-regulation)
* [Bcl6-Mediated Transcriptional Regulation of Follicular Helper T cells (TFH)](https://www.cell.com/trends/immunology/fulltext/S1471-4906(21)00026-0)
* [A distinct subpopulation of CD25− T-follicular regulatory cells localizes in the germinal centers](https://www.pnas.org/doi/abs/10.1073/pnas.1705551114)
* [T Follicular Helper Cell Biology: A Decade of Discovery and Diseases](https://www.cell.com/immunity/fulltext/S1074-7613(19)30191-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1074761319301918%3Fshowall%3Dtrue)
* [Innate Lymphoid Cells: 10 Years On](https://www.sciencedirect.com/science/article/pii/S0092867418309115)
* [Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors](https://www.science.org/doi/10.1126/science.aah4573?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
* [A cell atlas of human thymic development defines T cell repertoire formation](https://www.science.org/doi/10.1126/science.aay3224?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
* [Deciphering the fate of slan+-monocytes in human tonsils by gene expression profiling](https://faseb.onlinelibrary.wiley.com/doi/abs/10.1096/fj.202000181R)
* [FDC-SP, a Novel Secreted Protein Expressed by Follicular Dendritic Cells](https://www.jimmunol.org/content/169/5/2381.long)
















