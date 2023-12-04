#!/bin/bash

    
#SBATCH --job-name="q0qp6qyu_nkrug5p0"
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --time=10:00:00
#SBATCH -c 24
#SBATCH -n 1
#SBATCH -q normal
#SBATCH -p genD
#SBATCH --mem=120G
#SBATCH --error=./log/%x_%J.err
#SBATCH --output=./log/%x_%J.out


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/singlecell/software/cellranger/7.1.0/cellranger count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_162/jobs/q0qp6qyu_nkrug5p0/fastq --id q0qp6qyu_nkrug5p0 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A --include-introns=false;


echo [`date "+%Y-%m-%d %T"`] job finished
