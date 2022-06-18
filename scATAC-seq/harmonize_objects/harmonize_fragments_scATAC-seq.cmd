#!/bin/bash
#SBATCH --job-name="harmon_frags"
#SBATCH --output=./log/harmon_frags_%A_%a.out
#SBATCH --error=./log/harmon_frags_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --array=1-12

WD=/scratch/devel/rmassoni/tonsil_atlas/current
PATH_DATA=${WD}/multiome/1-cellranger_mapping/data/tonsil_atlas_metadata_multiome.csv
LINE=$(grep ATAC $PATH_DATA| cut -d',' -f1-2 | sed -n "$SLURM_ARRAY_TASK_ID"p)
subproject=$(echo $LINE | cut -d',' -f1)
gem_id=$(echo $LINE | cut -d',' -f2)
echo $subproject
echo $gem_id


bash harmonize_fragments_scATAC-seq.sh $subproject $gem_id