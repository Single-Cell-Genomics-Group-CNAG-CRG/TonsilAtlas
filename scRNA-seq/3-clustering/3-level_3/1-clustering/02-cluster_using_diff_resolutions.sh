for cell_type in NBC_MBC GCBC CD4_T Cytotoxic PC;
do
  echo $cell_type;
  sbatch --reservation=massoni-reserv --err "log/${cell_type}_clustering_diff_resolutions.err" --out "log/${cell_type}_clustering_diff_resolutions.log" -J "${cell_type}" -c 12 --time 23:00:00 --wrap "module purge; module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0; echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; Rscript ./02-cluster_using_diff_resolutions.R ${cell_type}; echo [`date "+%Y-%m-%d %T"`] job finished" 
done
