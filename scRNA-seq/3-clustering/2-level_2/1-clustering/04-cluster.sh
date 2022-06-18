PATH_OBJECTS="../../../results/R_objects/level_2/"
cell_types=$(ls $PATH_OBJECTS)
for cell_type in $cell_types;
do
    echo $cell_type;
    sbatch --reservation=massoni-reserv -J "${cell_type}" --err "./log/${cell_type}_clustering.err" --out "./log/${cell_type}_clustering.log" -c 8 --time 01:00:00 --wrap \
    "echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; \
    module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0 PANDOC; \
    R -e \"library(BiocStyle); \
          path_to_knit <- '/scratch/devel/rmassoni/tonsil_atlas/current/'; \
          cell_type <- '${cell_type}'; \
          rmarkdown::render( \
             '04-cluster.Rmd', \
             output_file = paste('reports/04-', cell_type, '_clustering_level_2.html', sep = ''), \
             knit_root_dir = path_to_knit \
          )\";
    echo [`date "+%Y-%m-%d %T"`] job finished"
done
