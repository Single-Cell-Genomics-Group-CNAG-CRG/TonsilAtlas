# for some reason it does not run when executed with bash script.sh. Copy and paste code to terminal
# Also, run "conda activate scrublet" before executing the code below
for subproject in ../../1-cellranger_mapping/projects/*; 
do 
	subproject=$(echo $subproject | rev | cut -d'/' -f1 | rev);
	for gem_id in ../../1-cellranger_mapping/projects/${subproject}/jobs/*;
	do 
		gem_id=$(echo $gem_id | rev | cut -d'/' -f1 | rev);
		sbatch -J "${gem_id}" --err "./log/${gem_id}.err" --out "./log/${gem_id}.log" -c 2 --time 00:15:00 --wrap \
			"echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; \
			python 01-doublet_detection_scrublet.py ${subproject} ${gem_id} ;
                        module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0 PANDOC;
			R -e \"library(BiocStyle); \
			      path_to_knit <- '/scratch/devel/rmassoni/tonsil_atlas/current/'; \
			      subproject <- '${subproject}'; \
			      gem_id <- '${gem_id}'; \
			      rmarkdown::render( \
			        '02-scrublet_results_summary.Rmd', \
			        output_file = paste('reports/scrublet_results_summary-', '${subproject}', '-', '${gem_id}', '.html', sep = ''), \
			        knit_root_dir = path_to_knit \
			      )\";
	      	echo [`date "+%Y-%m-%d %T"`] job finished"
		echo $subproject;
		echo $gem_id;
	done;
done
