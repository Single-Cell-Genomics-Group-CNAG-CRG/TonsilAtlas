# An Atlas of Cells in the Human Tonsil 

Palatine tonsils are secondary lymphoid organs (SLOs) representing the first line of immunological defense
against inhaled or ingested pathogens. We generated an atlas of the human tonsil composed of >556,000 cells
profiled across five different data modalities, including single-cell transcriptome, epigenome, proteome, and
immune repertoire sequencing, as well as spatial transcriptomics. This census identified 121 cell types and states,
defined developmental trajectories, and enabled an understanding of the functional units of the tonsil. Exemplarily,
we stratified myeloid slan-like subtypes, established a BCL6 enhancer as locally active in follicle-associated T
and B cells, and identified SIX5 as putative transcriptional regulator of plasma cell maturation. Analyses of a
validation cohort confirmed the presence, annotation, and markers of tonsillar cell types and provided evidence
of age-related compositional shifts. We demonstrate the value of this resource by annotating cells from B cell-derived
mantle cell lymphomas, linking transcriptional heterogeneity to normal B cell differentiation states of the human tonsil.

This repository contains all the scripts, notebooks and reports to reproduce all analysis from our manuscript
entitled [_An Atlas of Cells in the Human Tonsil_](https://www.cell.com/immunity/fulltext/S1074-7613(24)00031-1), published in Immunity in 2024.
Here, we describe how to access the data, document the most important packages and versions used, and explain how
to navigate the directories and files in this repository.

![](data/TonsilAtlasPic.png)


## Data

The data has been deposited in five levels of organization, from raw to processed data:


* Level 1: raw data. All fastq files for all data modalities have been deposited at ArrayExpress under accession id [E-MTAB-13687](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13687).
* Level 2: matrices. All data modalities correspond to different technologies from 10X Genomics. As such, they were mapped with different flavors of CellRanger (CR). The most important files in the ‘‘outs’’ folder of every CR run (including all matrices) have been deposited in [Zenodo](https://zenodo.org/records/10373041).
* Level 3: Seurat Objects. All data was analyzed within the Seurat ecosystem. We have archived in [Zenodo](https://zenodo.org/records/8373756) all Seurat Objects that contain the raw and processed counts, dimensionality reductions (PCA, Harmony, UMAP), and metadata needed to reproduce all figures from this manuscript.
* Level 4: to allow for programmatic and modular access to the whole tonsil atlas dataset, we developed [HCATonsilData](https://bioconductor.org/packages/release/data/experiment/html/HCATonsilData.html), available on BioConductor. HCATonsilData provides a vignette which documents how to navigate and understand the data. It also provides access to the glossary to traceback all annotations in the atlas. In addition, we will periodically update the annotations as we refine it with suggestions from the community.
* Level 5: interactive mode. Our tonsil atlas has been included as a reference in [Azimuth](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Tonsil%20v2), which allows interactive exploration of cell type markers on the web.

We refer to the READMEs in the Zenodo repositories for an explanation of how to access the matrices and Seurat objects. We have a separate repository ([TonsilAtlasCAP](https://github.com/massonix/TonsilAtlasCAP)) with scripts and documentation to download and remap all the fastq files from [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13687)


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
* 4-revision: we include one folder with the analysis performed to answer each of the major reviews during the revision of the manuscript.
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
















