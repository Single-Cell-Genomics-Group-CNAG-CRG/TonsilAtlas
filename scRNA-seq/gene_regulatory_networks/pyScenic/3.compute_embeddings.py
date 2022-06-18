import argparse
import os
import umap
from MulticoreTSNE import MulticoreTSNE as TSNE
import pandas as pd
import loompy as lp

################################
# Define arguments             #
################################
parser = argparse.ArgumentParser()
parser.add_argument("--PROJECT_NAME", type=str, help='name of the project that will be used for base folder and the output files')
parser.add_argument("--BASE_FOLDERNAME", type=str, help='base folder path where to save the object')
args = parser.parse_args()

project_name = args.PROJECT_NAME
base_foldername = args.BASE_FOLDERNAME


################################################################################
################################################################################
##   0. Preparing the parameters and files needed for the pySCENIC pipeline    #
################################################################################
################################################################################

################################
# Folder structure             #
################################
RESOURCES_FOLDERNAME = os.path.join(base_foldername, '{}/data/'.format(project_name))
RESULTS_FOLDERNAME = os.path.join(base_foldername, '{}/results/pyscenic_output/'.format(project_name))
FIGURES_FOLDERNAME = os.path.join(base_foldername, '{}/figures/'.format(project_name))

################################
# Input files                #
################################
FINAL_LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}_scenic_output.loom'.format(project_name))

################################
# Result files                #
################################
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}_auc.csv'.format(project_name))
REGULONS_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}_regulons.csv'.format(project_name))
AUCELL_TSNE = os.path.join(RESULTS_FOLDERNAME, '{}_aucell_tsne.tsv'.format(project_name))
AUCELL_UMAP = os.path.join(RESULTS_FOLDERNAME, '{}_aucell_umap.tsv'.format(project_name))

##################################################################################
##################################################################################
##                  3. NON-LINEAR PROJECTION AND CLUSTERING                      #
##################################################################################
##################################################################################

# how to work with loom files https://linnarssonlab.org/loompy/apiwalkthrough/index.html

lf = lp.connect(FINAL_LOOM_FNAME, mode='r+', validate=False )
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene)
lf.close()
auc_mtx.to_csv(AUCELL_MTX_FNAME, header=True)
regulons.to_csv(REGULONS_MTX_FNAME, header=True)


# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['AUCell_UMAP1', 'AUCell_UMAP2'], index=auc_mtx.index).to_csv(AUCELL_UMAP, sep='\t')

# tSNE
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['AUCell_tSNE1', 'AUCell_tSNE2'], index=auc_mtx.index).to_csv( AUCELL_TSNE, sep='\t')
