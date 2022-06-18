#!/bin/bash
    

#SBATCH --job-name=yja17lfh_kdslhq0g
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/yja17lfh_kdslhq0g_%x_%J.err
#SBATCH --output=./log/yja17lfh_kdslhq0g_%x_%J.out
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

/scratch/groups/hheyn/software/cellranger-atac/1.2.0/cellranger-atac count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scATAC_seq/1-cellranger_mapping/projects/BCLLATLAS_26/jobs/yja17lfh_kdslhq0g/fastq --id yja17lfh_kdslhq0g --sample yja17lfh_kdslhq0g --localcores 24 --localmem 64 --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-atac-GRCh38-1.2.0;

echo [`date "+%Y-%m-%d %T"`] job finished
