#!/bin/bash
    

#SBATCH --job-name=dhqzxpdk_h8em75r9
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/dhqzxpdk_h8em75r9_%x_%J.err
#SBATCH --output=./log/dhqzxpdk_h8em75r9_%x_%J.out
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=14


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

/scratch/groups/hheyn/software/cellranger-atac/1.2.0/cellranger-atac count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scATAC_seq/1-cellranger_mapping/projects/BCLLATLAS_30/jobs/dhqzxpdk_h8em75r9/fastq --id dhqzxpdk_h8em75r9 --sample dhqzxpdk_h8em75r9 --localcores 24 --localmem 64 --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-atac-GRCh38-1.2.0;

echo [`date "+%Y-%m-%d %T"`] job finished