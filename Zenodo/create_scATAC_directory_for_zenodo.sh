#!/bin/bash

# Define the paths and subprojects
base_path="/scratch/devel/rmassoni/tonsil_atlas/current/scATAC-seq/1-cellranger_mapping/projects"
destination_path="/scratch/devel/rmassoni/zenodo_20240525/scATAC-seq"
subprojects="BCLLATLAS_17 BCLLATLAS_20 BCLLATLAS_26 BCLLATLAS_30 BCLLATLAS_49"

# Function to compare MD5 checksums
compare_md5sum() {
    local src_file=$1
    local dst_file=$2

    local src_md5=$(md5sum "$src_file" | awk '{ print $1 }')
    local dst_md5=$(md5sum "$dst_file" | awk '{ print $1 }')

    if [ "$src_md5" == "$dst_md5" ]; then
        echo "MD5 checksum match for $dst_file"
    else
        echo "MD5 checksum mismatch for $dst_file"
    fi
}

# Change to the base path
cd $base_path

# Iterate through each subproject
for subproject in $subprojects; do
    echo "Processing subproject: $subproject"
    gem_ids=$(ls ${subproject}/jobs/)
    
    # Iterate through each gem_id within the subproject
    for gem_id in $gem_ids; do
        echo "Processing gem_id: $gem_id"
        
        # Define the source and target directories
        src_dir="${subproject}/jobs/${gem_id}/${gem_id}/outs"
        target_dir="${destination_path}/${subproject}/${gem_id}"
        
        # Create the target directory if it doesn't exist
        mkdir -p $target_dir
        
        # List of files to copy
        files_to_copy=(
            "web_summary.html"
            "summary.csv"
            "singlecell.csv"
            "peaks.bed"
            "peak_annotation.tsv"
            "fragments.tsv.gz.tbi"
            "fragments.tsv.gz"
            "filtered_peak_bc_matrix.h5"
        )
        
        # Copy each file to the target directory and verify MD5 checksum
        for file in "${files_to_copy[@]}"; do
            if [ -f "${src_dir}/${file}" ]; then
                cp "${src_dir}/${file}" "${target_dir}/"
                echo "Copied ${file} to ${target_dir}/"
                compare_md5sum "${src_dir}/${file}" "${target_dir}/${file}"
            else
                echo "File ${file} not found in ${src_dir}/"
            fi
        done
        
        # Copy the filtered_peak_bc_matrix directory recursively and verify MD5 checksum for each file
        if [ -d "${src_dir}/filtered_peak_bc_matrix" ]; then
            cp -r "${src_dir}/filtered_peak_bc_matrix" "${target_dir}/"
            echo "Copied filtered_peak_bc_matrix directory to ${target_dir}/"
            
            for file in $(find "${src_dir}/filtered_peak_bc_matrix" -type f); do
                rel_path=${file#${src_dir}/}
                compare_md5sum "$file" "${target_dir}/${rel_path}"
            done
        else
            echo "Directory filtered_peak_bc_matrix not found in ${src_dir}/"
        fi
    done
done
