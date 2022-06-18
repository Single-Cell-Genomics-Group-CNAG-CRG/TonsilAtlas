#!/bin/bash
    
#SBATCH --job-name="hg9au5oq_9ih1c0vd"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/hg9au5oq_9ih1c0vd_%x_%J.err
#SBATCH --output=./log/hg9au5oq_9ih1c0vd_%x_%J.out
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_25/jobs/hg9au5oq_9ih1c0vd/fastq --id hg9au5oq_9ih1c0vd --chemistry SC3Pv3 --expect-cells 5000 --localcores 24 --localmem 64 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;


echo [`date "+%Y-%m-%d %T"`] job finished