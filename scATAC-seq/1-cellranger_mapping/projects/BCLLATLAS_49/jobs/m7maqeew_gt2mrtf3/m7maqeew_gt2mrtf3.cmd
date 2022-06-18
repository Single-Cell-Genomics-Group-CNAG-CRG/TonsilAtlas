#!/bin/bash
    

#SBATCH --job-name=m7maqeew_gt2mrtf3
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/m7maqeew_gt2mrtf3_%x_%J.err
#SBATCH --output=./log/m7maqeew_gt2mrtf3_%x_%J.out
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=14


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

/scratch/groups/hheyn/software/cellranger-atac/1.2.0/cellranger-atac count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scATAC_seq/1-cellranger_mapping/projects/BCLLATLAS_49/jobs/m7maqeew_gt2mrtf3/fastq --id m7maqeew_gt2mrtf3 --sample m7maqeew_gt2mrtf3 --localcores 24 --localmem 64 --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-atac-GRCh38-1.2.0;

echo [`date "+%Y-%m-%d %T"`] job finished