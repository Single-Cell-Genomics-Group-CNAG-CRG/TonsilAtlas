# A periodic table of human tonsillar cells

Palatine tonsils are secondary lymphoid organs representing the first line of immunological defense against inhaled or ingested pathogens. Here, we provide a comprehensive census of cell types forming the human tonsil by applying single cell transcriptome, epigenome, proteome and adaptive immune repertoire as well as spatial transcriptomics, resulting in an atlas of >357,000 cells. We provide a glossary of 137 annotated cell types and states, and disentangle gene regulatory mechanisms that drive cells through specialized lineage trajectories. Exemplarily, we stratify multiple tonsil-resident myeloid slancyte subtypes, establish a distant BCL6 superenhancer as locally active in both follicle-associated B and T cells, and describe SIX5 as a potentially novel transcriptional regulator of plasma cell maturation. Further, our atlas acts as a reference map to pinpoint alterations observed in disease, shown here discovering immune-phenotype plasticity in tumoral cells and microenvironment shifts of mantle cell lymphomas (MCL). To facilitate such reference-based analysis, we develop HCATonsilData and SLOcatoR, a computational framework that provides programmatic and modular access to our dataset; and allows the straightforward annotation of future single-cell profiles from secondary lymphoid organs.


This repository contains all the scripts, notebooks and reports to reproduce all analysis from our recent preprint titled "A periodic table of tonsillar cells" (TODO: update title). Here, we describe how to access the data, document the most important packages and versions used, and explain how to navigate the directories and files in this repository.


## Data


The data has been deposited in five levels of organization, from raw to processed data:


* Level 1: raw data. All fastq files for all data modalities have been deposited at the European Genome Archive (EGA) under accession id XXX.
* Level 2: matrices. All data modalities correspond to different technologies from 10X Genomics. As such, they were mapped with different flavors of CellRanger (CR). The most important files in the “outs” folder of every CR run (including all matrices) have been deposited in Zenodo (TODO: include link).
* Level 3: Seurat Objects. All data was analyzed within the Seurat ecosystem. We have archived in [Zenodo](https://doi.org/10.5281/zenodo.6340174) all Seurat Objects that contain the raw and processed counts, dimensionality reductions (PCA, Harmony, UMAP), and metadata needed to reproduce all figures from this manuscript.
* Level 4: to allow for programmatic and modular access to the whole tonsil atlas dataset, we developed HCATonsilData, available on [GitHub](https://github.com/massonix/HCATonsilData). HCATonsilData provides a vignette which documents how to navigate and understand the data. In addition, we will periodically update the annotations as we refine it with suggestions from the  community.
* Level 5: interactive mode. We provide different iSEE instances to browse the data interactively: [http://shiny.imbei.uni-mainz.de:3838/iSEE_TonsilDataAtlas/](http://shiny.imbei.uni-mainz.de:3838/iSEE_TonsilDataAtlas/).


We refer to the READMEs in this repositories for a full explanation of the dataset.



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
* [chromVar](https://bioconductor.org/packages/release/bioc/html/chromVAR.html)
* [JASPAR2020](http://bioconductor.org/packages/release/data/annotation/html/JASPAR2020.html)
* [chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs)
* [Vireo 0.5.0](https://github.com/single-cell-genetics/vireo)
* [Scirpy 0.7.0](https://scverse.org/scirpy/latest/tutorials.html)
* [SPOTlight](https://bioconductor.org/packages/release/bioc/html/SPOTlight.html)
* [MAGIC](https://github.com/KrishnaswamyLab/MAGIC)
* [SPATA2](https://github.com/theMILOlab/SPATA2)



You can check the versions of other packages at the "Session Information" section of each html report.


# File system and name scheme

TODO


# Relevant literature

TODO
















