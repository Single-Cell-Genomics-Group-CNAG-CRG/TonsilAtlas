#!/bin/bash
now=$(date)
echo $now
echo "Running pySCENIC script"

################################
# ACTIVATE CONDA             #
################################

# Activate conda environment with the right dependencies
## Get base path
CONDA_BASE=$(conda info --base)
## reference your conda installation path
source $CONDA_BASE/etc/profile.d/conda.sh
## Activate conda
conda activate scenic_protocol
export HDF5_USE_FILE_LOCKING='FALSE'

################################
# CREATE THE FOLDERS             #
################################

project_name="name_of_the_project"
base_folder="{path_to_project_folder}"

echo $project_name

################################
# Folder structure             #
################################
RESOURCES_FOLDERNAME=$base_folder$project_name"/data/"
RESULTS_FOLDERNAME=$base_folder$project_name"/results/pyscenic_output/"

AUXILLIARIES_FOLDERNAME="{path_to_auxilliaries_folder}" #folder with Auxilliary files (TF, rankingDB and motif annotations)

################################
# Input files                  #
################################
INITIAL_LOOM_SEURATOBJECT=$RESOURCES_FOLDERNAME$project_name".loom"


################################
# Files created                #
################################
ADJACENCIES_FNAME=$RESULTS_FOLDERNAME$project_name"_adjacencies.tsv"
MOTIFS_FNAME=$RESULTS_FOLDERNAME$project_name"_motifs.csv"
FINAL_LOOM_FNAME=$RESULTS_FOLDERNAME$project_name"_scenic_output.loom"

################################
# Auxilliary files             #
################################
# Downloaded fromm pySCENIC github repo: https://github.com/aertslab/pySCENIC/tree/master/resources
HUMAN_TFS_FNAME=$AUXILLIARIES_FOLDERNAME"TF/allTFs_hg38.txt"
# Ranking databases. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
RANKING_DBS_FNAMES=$AUXILLIARIES_FOLDERNAME"RankingDB/hg38/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
# Motif annotations. Downloaded from cisTargetDB: https://resources.aertslab.org/cistarget/
MOTIF_ANNOTATIONS_FNAME=$AUXILLIARIES_FOLDERNAME"Annotation/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"


################################################################################
################################################################################
##                           2. SCENIC STEPS                                   #
################################################################################
################################################################################

########################################################################################################
########################################################################################################
## STEP 1: Gene regulatory network inference, and generation of co-expression modules  ● Timing 95 min #
########################################################################################################
########################################################################################################

####
# GRN inference using the GRNBoost2 algorithm
####
#
echo "Running GRN step"
arboreto_with_multiprocessing.py \
  $INITIAL_LOOM_SEURATOBJECT \
  $HUMAN_TFS_FNAME \
  --method grnboost2 \
  --output $ADJACENCIES_FNAME \
  --num_workers 23 \
  --seed 777

##################################################################################
##################################################################################
## STEP 2-3: Candidate regulon generation and regulon prediction ● Timing 8 min  #
##################################################################################
##################################################################################

###################################
# Regulon generation & prediction #
###################################
echo "Running CTX step"
pyscenic ctx $ADJACENCIES_FNAME $RANKING_DBS_FNAMES \
  --annotations_fname $MOTIF_ANNOTATIONS_FNAME \
  --expression_mtx_fname $INITIAL_LOOM_SEURATOBJECT \
  --output $MOTIFS_FNAME \
  --num_workers 23

#########################################################################
#########################################################################
## STEP 4: Cellular enrichment (aka AUCell) from CLI ● Timing 1.7 min   #
#########################################################################
#########################################################################

################################
# Cellular enrichment (AUCELL) #
################################
echo "Running AUCell step"
pyscenic aucell \
    $INITIAL_LOOM_SEURATOBJECT \
    $MOTIFS_FNAME \
    --output $FINAL_LOOM_FNAME \
    --num_workers 23

#########################################################################
#########################################################################
## STEP 5: Compute Embedding   #
#########################################################################
#########################################################################
python3.6 3.compute_embeddings.py --PROJECT_NAME $project_name --BASE_FOLDERNAME $base_folder


#########################################################################
#########################################################################
## STEP 6: Compute RSS   #
#########################################################################
#########################################################################
variable_rss="annotation_level_1" #column from the metadata from which we want to compute the RSS
python3.6 4.compute_rss.py --PROJECT_NAME $project_name --VARIABLE_TO_RSS $variable_rss --BASE_FOLDERNAME $base_folder


now=$(date)
echo $now
