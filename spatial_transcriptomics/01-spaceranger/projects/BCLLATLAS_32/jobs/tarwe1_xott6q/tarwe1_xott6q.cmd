#!/bin/bash
    
    

module load PYTHON/2.7.5 
module load lims/1.2 

/scratch/devel/melosua/phd/10x_software/spaceranger/1.1.0/spaceranger count --id=tarwe1_xott6q          --transcriptome=/scratch/devel/melosua/phd/10x_software/refdata-gex-GRCh38-2020-A/          --fastqs=../projects/BCLLATLAS_32/fastq          --sample=tarwe1_xott6q          --image=../img/tarwe1_xott6q_V19S23-039_A1.jpg          --slide=V19S23-039          --area=A1          --localcores=8          --localmem=64          --slidefile=../data/spot_layout/V19S23-039.gpr

