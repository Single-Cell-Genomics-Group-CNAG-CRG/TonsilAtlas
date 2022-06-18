PATH_OBJECTS="../results/R_objects/level_2/"
cell_types=$(ls $PATH_OBJECTS)
for cell_type in $cell_types;
do
  echo $cell_type;
 sbatch --reservation=massoni-reserv -J "${cell_type}" --err "./log/${cell_type}_compositional_analysis.err" --out "./log/${cell_type}_compositional_analysis.log" -c 8 --time 01:00:00 --wrap \
 "echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; \
 module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0 PANDOC; \
 Rscript 01-differential_expression_analysis.R ${cell_type}
 echo [`date "+%Y-%m-%d %T"`] job finished"
done
