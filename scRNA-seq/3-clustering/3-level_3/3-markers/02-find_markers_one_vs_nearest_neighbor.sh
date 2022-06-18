PATH_OBJECTS="../../../results/R_objects/level_2/"
cell_types=$(ls $PATH_OBJECTS)
for cell_type in $cell_types;
do
  echo $cell_type;
    sbatch --reservation=massoni-reserv --err "log/markers_${cell_type}_one_vs_nn.err" --out "log/markers_${cell_type}_one_vs_nn.log" -J "${cell_type}" -c 24 --time 23:00:00 --wrap "module purge; module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0; echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; Rscript ./02-find_markers_one_vs_nearest_neighbor.R ${cell_type}; echo [`date "+%Y-%m-%d %T"`] job finished" 
done

