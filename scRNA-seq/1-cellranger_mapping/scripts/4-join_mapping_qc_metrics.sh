#!/bin/bash

####### non-hashed libraries #######

# generate paths to metrics.csv
outs_list=$(awk -F, '/not_hashed/ {print $1,$2}' data/tonsil_atlas_metadata.csv | \
    sort | \
    uniq | \
    awk 'BEGIN {OFS = "/"}; {print "projects", $1, "jobs", $2, $2, "outs"}')

# write joint file
header=$(for outs_dir in $outs_list; do
    if [[ -d $outs_dir ]]; then
        tmp=$(head -n1 "${outs_dir}/metrics_summary.csv")
        echo "gem_id,${tmp}"
        break
    fi
done) 
echo $header > ../results/tables/cellranger_mapping/cellranger_mapping_metrics_not_hashed.csv

for outs_dir in $outs_list; do
    if [[ -d $outs_dir ]]; then
        gem_id=$(echo $outs_dir | awk -F/ '{print $4}')
        metrics=$(egrep -v "Estimated Number of Cells" "${outs_dir}/metrics_summary.csv")
        echo "$gem_id,$metrics" >> ../results/tables/cellranger_mapping/cellranger_mapping_metrics_not_hashed.csv
    fi
done


####### hashed libraries #######

# generate paths to metrics.csv
outs_list=$(awk -F, '!/not_hashed/ {print $1,$2}' data/tonsil_atlas_metadata.csv | \
    sort | \
    uniq | \
    awk 'BEGIN {OFS = "/"}; {print "projects", $1, "jobs", $2, $2, "outs"}')

# write joint file
header=$(for outs_dir in $outs_list; do
    if [[ -d $outs_dir ]]; then
        tmp=$(head -n1 "${outs_dir}/metrics_summary.csv")
        echo "gem_id,${tmp}"
        break
    fi
done)
echo $header > ../results/tables/cellranger_mapping/cellranger_mapping_metrics_hashed.csv

for outs_dir in $outs_list; do
    if [[ -d $outs_dir ]]; then
        gem_id=$(echo $outs_dir | awk -F/ '{print $4}')
        metrics=$(egrep -v "Estimated Number of Cells" "${outs_dir}/metrics_summary.csv")
        echo "$gem_id,$metrics" >> ../results/tables/cellranger_mapping/cellranger_mapping_metrics_hashed.csv
    fi
done
