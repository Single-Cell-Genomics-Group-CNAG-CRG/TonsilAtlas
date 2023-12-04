#!/bin/bash
    
#SBATCH --job-name="rssxkrxg_edqhwe1o"
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/rssxkrxg_edqhwe1o_%x_%J.err
#SBATCH --output=./log/rssxkrxg_edqhwe1o_%x_%J.out
#SBATCH --time=10:00:00
#SBATCH -c 24
#SBATCH -n 1
#SBATCH -q normal
#SBATCH -p genD
#SBATCH --mem=120G


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/singlecell/software/cellranger/7.1.0/cellranger count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_162/jobs/rssxkrxg_edqhwe1o/fastq --id rssxkrxg_edqhwe1o --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A --include-introns=false;


echo [`date "+%Y-%m-%d %T"`] job finished
