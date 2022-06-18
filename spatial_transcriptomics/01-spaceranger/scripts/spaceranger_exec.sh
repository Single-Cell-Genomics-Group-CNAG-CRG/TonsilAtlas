#!/bin/bash

# module load python3/3.6.5
# Activate conda environment with the right dependencies
## Get base path
CONDA_BASE=$(conda info --base)
## reference your conda installation path
source $CONDA_BASE/etc/profile.d/conda.sh
## Activate conda
conda activate spaceranger

for projects in BCLLATLAS_51 # BCLLATLAS_32 
do
# Set common parameters
  pipe_fn="${projects}_pipe"
  # Create pipefile to append pipe commands
  rm -r $pipe_fn
  touch $pipe_fn
done

# Copy the fastq symbolic links to the project
python3 2-copy_lims_files_spaceranger.py info.txt ../projects ../data/sample_id.txt

declare -A subproj_dict=( \
  ["tarwe1_xott6q"]="BCLLATLAS_32" \
  ["c28w2r_7jne4i"]="BCLLATLAS_32" \
  ["esvq52_nluss5"]="BCLLATLAS_32" \
  ["p7hv1g_tjgmyj"]="BCLLATLAS_32" \
  ["gcyl7c_cec61b"]="BCLLATLAS_51" \
  ["zrt7gl_lhyyar"]="BCLLATLAS_51" \
  ["qvwc8t_2vsr67"]="BCLLATLAS_51" \
  ["exvyh1_66caqq"]="BCLLATLAS_51")

# c28w2r_7jne4i esvq52_nluss5 p7hv1g_tjgmyj
for id in gcyl7c_cec61b zrt7gl_lhyyar qvwc8t_2vsr67 exvyh1_66caqq # tarwe1_xott6q c28w2r_7jne4i esvq52_nluss5 p7hv1g_tjgmyj 
do
  subproj=${subproj_dict[${id}]}

  # printf "python 3-make_spaceranger.py --spaceranger /scratch/devel/melosua/phd/10x_software/spaceranger-1.1.0 --subproject {$subproj} --reference_path /scratch/devel/melosua/phd/10x_software/refdata-gex-GRCh38-2020-A/ --gem_id {$id} --metadata ../data/sample_id.txt\n"
  
  # 3-make_spaceranger.py writes a bash script ready to launch a job to the cluster
  python3 3-make_spaceranger.py \
    --spaceranger /scratch/groups/hheyn/software/spaceranger/1.1.0 \
    --subproject $subproj \
    --reference_path /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A/ \
    --gem_id $id \
    --metadata ../data/sample_id.txt \
  
  # Give execution permission to the file  
  chmod 755 ../projects/$subproj/jobs/$id/$id.cmd
  
  # Send job to the cluster
  ## Normal
  # printf "$id.spaceranger\t$id\t-\t.\tn\t08:00:00\t1\t10\t.\tmodule purge; ./../projects/$subproj/jobs/$id/$id.cmd\n" >> "$pipe_fn"
  ## Cluster mode
  printf "$id.spaceranger\t$id\t-\t.\tn\t04:00:00\t1\t1\t.\tmodule purge; ./../projects/$subproj/jobs/$id/$id.cmd\n" >> "$pipe_fn"
    
done

for projects in BCLLATLAS_51 # BCLLATLAS_32
do
  # Set common parameters
  pipe_fn="${projects}_pipe"
  /home/devel/melosua/bin/cnag_pipeline.pl "$pipe_fn"
done
