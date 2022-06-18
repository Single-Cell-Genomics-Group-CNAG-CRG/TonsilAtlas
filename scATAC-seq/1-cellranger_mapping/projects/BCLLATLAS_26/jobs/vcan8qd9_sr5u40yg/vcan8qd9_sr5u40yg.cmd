#!/bin/bash
    

#SBATCH --job-name=vcan8qd9_sr5u40yg
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/vcan8qd9_sr5u40yg_%x_%J.err
#SBATCH --output=./log/vcan8qd9_sr5u40yg_%x_%J.out
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=14


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

/scratch/groups/hheyn/software/cellranger-atac/1.2.0/cellranger-atac count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scATAC_seq/1-cellranger_mapping/projects/BCLLATLAS_26/jobs/vcan8qd9_sr5u40yg/fastq --id vcan8qd9_sr5u40yg --sample vcan8qd9_sr5u40yg --localcores 24 --localmem 64 --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-atac-GRCh38-1.2.0;

echo [`date "+%Y-%m-%d %T"`] job finished