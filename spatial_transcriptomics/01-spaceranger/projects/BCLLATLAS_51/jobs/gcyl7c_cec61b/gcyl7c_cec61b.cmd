#!/bin/bash
    
    

module load PYTHON/2.7.5 
module load lims/1.2 

/scratch/groups/hheyn/software/spaceranger/1.1.0//spaceranger count --id=gcyl7c_cec61b          --transcriptome=/scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A/          --fastqs=../projects/BCLLATLAS_51/fastq          --sample=gcyl7c_cec61b          --image=../img/gcyl7c_cec61b_V10M16-059_A1.jpg          --slide=V10M16-059          --area=A1          --localcores=2          --slidefile=../data/spot_layout/V10M16-059.gpr 	 --jobmode=/scratch/groups/hheyn/software/spaceranger/1.1.0//external/martian/jobmanagers/slurm.template

