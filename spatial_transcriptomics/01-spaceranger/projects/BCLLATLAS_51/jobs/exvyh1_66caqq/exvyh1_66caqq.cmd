#!/bin/bash
    
    

module load PYTHON/2.7.5 
module load lims/1.2 

/scratch/groups/hheyn/software/spaceranger/1.1.0//spaceranger count --id=exvyh1_66caqq          --transcriptome=/scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A/          --fastqs=../projects/BCLLATLAS_51/fastq          --sample=exvyh1_66caqq          --image=../img/exvyh1_66caqq_V10M16-059_D1.jpg          --slide=V10M16-059          --area=D1          --localcores=2          --slidefile=../data/spot_layout/V10M16-059.gpr 	 --jobmode=/scratch/groups/hheyn/software/spaceranger/1.1.0//external/martian/jobmanagers/slurm.template

