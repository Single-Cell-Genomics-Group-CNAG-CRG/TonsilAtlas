for i in $(grep hashed_cdna ../../1-cellranger_mapping/data/tonsil_atlas_metadata.csv);
do
	subproject=$(echo $i | cut -d',' -f1);
	gem_id=$(echo $i | cut -d',' -f2);
	echo $subproject;
	echo $gem_id;
	sbatch -J "${gem_id}" --err "./log/${gem_id}.err" --out "./log/${gem_id}.log" -c 3 --time 00:35:00 --wrap \
	"echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; \
    module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0 PANDOC; \
	R -e \"library(BiocStyle); \
	      path_to_knit <- '/scratch/devel/rmassoni/tonsil_atlas/current/'; \
	      subproject <- '${subproject}'; \
	      gem_id <- '${gem_id}'; \
	      save_object_path <- 'scRNA-seq/results/R_objects/demultiplexed/seurat_${gem_id}_demultiplexed.rds'; \
	      rmarkdown::render( \
	        '01-demultiplex_hashtags.Rmd', \
	        output_file = paste('reports/demultiplex_hashtags-', gem_id, '.html', sep = ''), \
	        knit_root_dir = path_to_knit \
	      )\";
  	echo [`date "+%Y-%m-%d %T"`] job finished"
done
