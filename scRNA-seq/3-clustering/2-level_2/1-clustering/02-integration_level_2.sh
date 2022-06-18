PATH_OBJECTS="../../../results/R_objects/level_2/"
cell_types=$(ls $PATH_OBJECTS)
for cell_type in $cell_types;
do
  echo $cell_type;
  sbatch --reservation=massoni-reserv --err "log/test_${cell_type}.err" --out "log/test_${cell_type}.log" -J "rf_${res}" -c 12 --time 23:00:00 --wrap "module purge; module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0; echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; Rscript ./02-integration_level_2.R ${cell_type}; echo [`date "+%Y-%m-%d %T"`] job finished" 
done
