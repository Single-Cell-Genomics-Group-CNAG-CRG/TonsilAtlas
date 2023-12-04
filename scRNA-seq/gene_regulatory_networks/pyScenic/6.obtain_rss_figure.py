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
project_name="name_of_the_project"
base_foldername="{path_to_project_folder}"
variable_to_rss="annotation_level_1"


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
AUCELL_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}_auc.csv'.format(project_name))

################################
# Result files                  #
################################
RSS_MTX = os.path.join(RESULTS_FOLDERNAME, '{}_{}_rss_matrix.csv'.format(project_name, variable_to_rss))


################################
##########################################################################################################################################################################################################################################
############################################################################################################################################################
################################################################################################################################################################################################################

def plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[1]
    data = rss.T[cell_type].sort_values(ascending=False)[0:max_n]
    ax.plot(np.arange(len(data)), data, '.')
    ax.set_ylim([np.floor(data.min()*100.0)/100.0, np.ceil(data.max()*100.0)/100.0])
    ax.set_ylabel('RSS')
    ax.set_xlabel('Regulon')
    ax.set_title(cell_type)
    ax.set_xticklabels([])
    font = {
        'color':  'red',
        'weight': 'normal',
        'size': 8,
    }
    for idx, (regulon_name, rss_val) in enumerate(zip(data[0:top_n].index, data[0:top_n].values)):
        ax.plot([idx, idx], [rss_val, rss_val], 'r.')
        ax.text(idx+(max_n/25), rss_val, regulon_name, fontdict=font, horizontalalignment='left', verticalalignment='center')


def savepdf(fname: str, fig, folder: str=FIGURES_FOLDERNAME) -> None:
    """
    Save figure as vector-based PDF image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname), format='pdf')




##################################################################################
##################################################################################
##            1. READ INPUT FILES               #
##################################################################################
##################################################################################

rss = pd.read_csv(RSS_MTX, index_col=0)

##################################################################################
##################################################################################
##            2. OBTAIN FIGURES               #
##################################################################################
##################################################################################
sns.set()
sns.set(style='whitegrid', font_scale=0.8)
fig, (ax1) = plt.subplots(1, 1, figsize=(5, 15), dpi=100)
plot_rss(rss, 'PC', top_n=20, ax=ax1)
plt.tight_layout()

##########################
# RSS for each celltype  #
##########################
savepdf('PC - rss for each celltype.pdf', fig)
