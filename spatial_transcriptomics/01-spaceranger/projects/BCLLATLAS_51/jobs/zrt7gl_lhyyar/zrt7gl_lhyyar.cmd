#!/bin/bash
    
    

module load PYTHON/2.7.5 
module load lims/1.2 

/scratch/groups/hheyn/software/spaceranger/1.1.0//spaceranger count --id=zrt7gl_lhyyar          --transcriptome=/scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A/          --fastqs=../projects/BCLLATLAS_51/fastq          --sample=zrt7gl_lhyyar          --image=../img/zrt7gl_lhyyar_V10M16-059_B1.jpg          --slide=V10M16-059          --area=B1          --localcores=2          --slidefile=../data/spot_layout/V10M16-059.gpr 	 --jobmode=/scratch/groups/hheyn/software/spaceranger/1.1.0//external/martian/jobmanagers/slurm.template

