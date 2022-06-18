#!/bin/bash


# Load modules
module purge
module load gcc/4.9.3-gold zlib/1.2.8 HTSLIB/latest


# Define variables
subproject=$1
gem_id=$2
WD=/scratch/devel/rmassoni/tonsil_atlas/current
OUTDIR=${WD}/data/raw_data_figures/scATAC_CD4_T/fragments
OUTFILE=${OUTDIR}/${gem_id}_atac_fragments_renamed_cells.tsv
echo $subproject;
echo $gem_id;
if [[ "$subproject" == "BCLLATLAS_43" ]];
then
    echo "experiment 1"
    frag_file=${WD}/multiome/1-cellranger_mapping/projects/experiment_1/jobs/${gem_id}/${gem_id}/outs/atac_fragments.tsv.gz
elif [[ "$subproject" == "BCLLATLAS_48" ]];
then
    echo "experiment 2"
    frag_file=${WD}/multiome/1-cellranger_mapping/projects/experiment_2/jobs/${gem_id}/${gem_id}/outs/atac_fragments.tsv.gz
fi


# Rename fragments files and reindex
zcat ${frag_file} | awk -v frag="$gem_id" 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,frag"_"$4,$5}' >| $OUTFILE
bgzip $OUTFILE
tabix -p bed "${OUTFILE}.gz"