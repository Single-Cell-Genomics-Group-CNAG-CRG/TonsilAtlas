#!/bin/bash

####### CITE-seq libraries #######

# generate paths to metrics.csv
outs_list=$(sed '1d' data/tonsil_atlas_metadata.csv | awk -F, '{print $1,$2}' | sort | uniq | awk 'BEGIN {OFS = "/"}; {print "projects", $1, "jobs", $2, $2, "outs","per_sample_outs" , $2}')

# write joint file
cnt=0
for outs_dir in $outs_list; do        
    if [[ -d $outs_dir ]]; then
        if [[ $cnt == 0 ]]; then
            head -1 "${outs_dir}/metrics_summary.csv" | awk -F"," '{print "subproject"",""gem"","$0}' > results/cellranger_mapping_metrics_cite_seq.csv
        fi
        gem_id=$(echo $outs_dir | awk -F/ '{print $4}')
        subproject=$(echo $outs_dir | awk -F/ '{print $2}');
        sed '1d' "${outs_dir}/metrics_summary.csv" | awk -v gem=$gem_id -v subp=$subproject '{print subp","gem","$0}' >> results/cellranger_mapping_metrics_cite_seq.csv
        cnt=$((cnt+1))          
    fi
done
