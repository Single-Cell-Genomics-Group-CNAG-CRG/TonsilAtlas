import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import zlib
import base64

import os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

from IPython.display import HTML, display

################################
# Define arguments             #
################################
parser = argparse.ArgumentParser()
parser.add_argument("--PROJECT_NAME", type=str, help='name of the project that will be used for base folder and the output files')
parser.add_argument("--VARIABLE_TO_RSS", type=str, help='variable (column name from metadata) to compute the rss')
parser.add_argument("--BASE_FOLDERNAME", type=str, help='base folder path where to save the object')
args = parser.parse_args()

project_name = args.PROJECT_NAME
variable_to_rss = args.VARIABLE_TO_RSS
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
FIGURES_FOLDERNAME = os.path.join(base_foldername, '{}/results/figures/rss_results/'.format(project_name))

################################
# Input files                  #
################################
METADATA_FNAME = os.path.join(RESOURCES_FOLDERNAME,  '{}_metadata.csv'.format(project_name))
#METADATA_FNAME = os.path.join(RESOURCES_FOLDERNAME,  'PC/metadata.csv'.format(project_name))
AUCELL_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}_auc.csv'.format(project_name))

################################
# Result files                  #
################################
RSS_MTX = os.path.join(RESULTS_FOLDERNAME, '{}_{}_rss_matrix.csv'.format(project_name, variable_to_rss))


################################


##################################################################################
##################################################################################
##            1. READ INPUT FILES               #
##################################################################################
##################################################################################

##########################
# AUC_matrix  #
##########################
#lf = lp.connect(FINAL_LOOM_FNAME, mode='r+', validate=False )
#auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
#lf.close()
auc_mtx = pd.read_csv(AUCELL_FNAME, index_col=0)

##########################
# METADATA  #
##########################
df_metadata = pd.read_csv(METADATA_FNAME, index_col = 0, encoding="latin-1")

# Specify the variable for defining the groups for the specificity score
groups_specificity_object = df_metadata[variable_to_rss]


##################################################################################
##################################################################################
##            3. COMPUTE RSS FOR CLUSTERS              #
##################################################################################
##################################################################################

#######################################################################
# ENSURE BOTH OBJECTS HAVE THE SAME ORDER OF CELLS #
#######################################################################
df = auc_mtx
index = df.index
a_list = list(index)
df["CellID"] = a_list

index = df_metadata.index
a_list = list(index)

# Create the dictionary that defines the order for sorting
sorterIndex = dict(zip(a_list, range(len(a_list))))


# Generate a rank column that will be used to sort
# the dataframe numerically
df['CellID_Rank'] = df['CellID'].map(sorterIndex)

# Here is the result asked with the lexicographic sort
# Result may be hard to analyze, so a second sorting is
# proposed next
## NOTE:
## Newer versions of pandas use 'sort_values' instead of 'sort'
df.sort_values(['CellID_Rank'],
        ascending = [True], inplace = True)
df.drop('CellID_Rank', 1, inplace = True)
df.drop('CellID', 1, inplace = True)
print(df)


#######################################################################
#                            COMPUTE RSS                              #
#######################################################################
rss = regulon_specificity_scores(df, groups_specificity_object)
rss.head()

#######################################################################
#                        WRITE IT INTO CSV FILE FOR R                 #
#######################################################################
rss.to_csv(RSS_MTX, header=True)
