#!/bin/bash
    
#SBATCH --job-name="jb6vuao4_g4vi9ur0"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/jb6vuao4_g4vi9ur0_%x_%J.err
#SBATCH --output=./log/jb6vuao4_g4vi9ur0_%x_%J.out
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_05/jobs/jb6vuao4_g4vi9ur0/fastq --id jb6vuao4_g4vi9ur0 --chemistry SC3Pv3 --expect-cells 5000 --localcores 24 --localmem 64 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;


echo [`date "+%Y-%m-%d %T"`] job finished