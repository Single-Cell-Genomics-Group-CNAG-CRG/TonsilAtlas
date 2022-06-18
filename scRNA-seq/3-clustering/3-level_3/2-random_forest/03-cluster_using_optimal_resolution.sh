for cell_type in NBC_MBC GCBC CD4_T Cytotoxic PC;
do
    echo $cell_type;
    sbatch --reservation=massoni-reserv -J "${cell_type}" --err "./log/${cell_type}_cluster.err" --out "./log/${cell_type}_cluster.log" -c 12 --time 01:00:00 --wrap \
    "echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; \
    module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0 PANDOC; \
    R -e \"library(BiocStyle); \
          path_to_knit <- '/scratch/devel/rmassoni/tonsil_atlas/current/'; \
          cell_type <- '${cell_type}'; \
          rmarkdown::render( \
             '03-cluster_using_optimal_resolution.Rmd', \
             output_file = paste('reports/03-', cell_type, '_cluster_using_optimal_resolution.html', sep = ''), \
             knit_root_dir = path_to_knit \
          )\";
    echo [`date "+%Y-%m-%d %T"`] job finished"
done
