#!/bin/bash
    
#SBATCH --job-name="izi9unx1_8qdzhivu"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/velocyto_izi9unx1_8qdzhivu_%J.err
#SBATCH --output=./log/velocyto_izi9unx1_8qdzhivu_%J.out
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

export HDF5_USE_FILE_LOCKING="FALSE"

REPEAT_MASK_GTF="/scratch/groups/hheyn/data/reference/GRCh38_repeat_mask/GRCh38_repeat_mask.gtf"
CELLRANGER_OUTPUT_DIR="/scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_05/jobs/izi9unx1_8qdzhivu/izi9unx1_8qdzhivu"
REFERENCE_GTF="/scratch/devel/rmassoni/reference/human/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"

velocyto run10x -m $REPEAT_MASK_GTF $CELLRANGER_OUTPUT_DIR $REFERENCE_GTF

echo [`date "+%Y-%m-%d %T"`] job finished