#!/bin/bash
    

#SBATCH --job-name=jbp85ve3_3179tpob
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/jbp85ve3_3179tpob_%x_%J.err
#SBATCH --output=./log/jbp85ve3_3179tpob_%x_%J.out
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

/scratch/groups/hheyn/software/cellranger-atac/1.2.0/cellranger-atac count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scATAC_seq/1-cellranger_mapping/projects/BCLLATLAS_17/jobs/jbp85ve3_3179tpob/fastq --id jbp85ve3_3179tpob --sample jbp85ve3_3179tpob --localcores 24 --localmem 64 --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-atac-GRCh38-1.2.0;

echo [`date "+%Y-%m-%d %T"`] job finished
