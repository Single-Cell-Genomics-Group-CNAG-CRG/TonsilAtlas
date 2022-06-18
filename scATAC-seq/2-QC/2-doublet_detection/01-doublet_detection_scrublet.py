#!/usr/bin/env python
# coding: utf-8

#####################################################################################
########################## Doublet Detection with Scrublet ##########################
#####################################################################################
 
# - author: "Paula Soler-Vila"
# - date: "12/02/2020"

# # Objective
# 
# This modified script reads the filtered_feature_bc_matrix from a specific GEM-id and runs scrublet.
# It saves a dataframe with the doublet score and prediction for each cell, the UMAPs with the
# scores projected and the histograms of their distribution.
# 

# Import packages
import numpy as np
import pandas as pd
import scipy.io
import sys
import os
import scrublet as scr
import pynndescent
import datetime
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

starting_time = datetime.datetime.now()
print("Starting job at " + starting_time.strftime("%Y-%m-%d %H:%M:%S"))

# Initialize variables
gem_id = sys.argv[1]


# These functions haven been modified from the https://github.com/AllonKleinLab/scrublet.

def custom_cmap(rgb_list):
    import matplotlib.pyplot as plt
    rgb_list = np.array(rgb_list)
    cmap = plt.cm.Reds
    cmap = cmap.from_list(rgb_list.shape[0],rgb_list)
    return cmap

def darken_cmap(cmap, scale_factor):
    cdat = np.zeros((cmap.N, 4))
    for ii in range(cdat.shape[0]):
        curcol = cmap(ii)
        cdat[ii,0] = curcol[0] * scale_factor
        cdat[ii,1] = curcol[1] * scale_factor
        cdat[ii,2] = curcol[2] * scale_factor
        cdat[ii,3] = 1
    cmap = cmap.from_list(cmap.N, cdat)
    return cmap

def plot_embedding(self, embedding_name, score='raw', marker_size=5, 
                   order_points=False, fig_size=(8,4), color_map=None):
    ''' Plot doublet predictions on 2-D embedding of observed transcriptomes '''

    #from matplotlib.lines import Line2D
    if embedding_name not in self._embeddings:
        print('Cannot find "{}" in embeddings. First add the embedding using `set_embedding`.'.format(embedding_name))
        return

    # TO DO: check if self.predicted_doublets exists; plot raw scores only if it doesn't

    fig, axs = plt.subplots(1, 2, figsize = fig_size)

    x = self._embeddings[embedding_name][:,0]
    y = self._embeddings[embedding_name][:,1]
    xl = (x.min() - x.ptp() * .05, x.max() + x.ptp() * 0.05)
    yl = (y.min() - y.ptp() * .05, y.max() + y.ptp() * 0.05)

    ax = axs[1]
    if score == 'raw':
        color_dat = self.doublet_scores_obs_
        vmin = color_dat.min()
        vmax = color_dat.max()
        if color_map is None:
            cmap_use = darken_cmap(plt.cm.Reds, 0.9)
        else:
            cmap_use = color_map
    elif score == 'zscore':
        color_dat = self.z_scores_
        vmin = -color_dat.max()
        vmax = color_dat.max()
        if color_map is None:
            cmap_use = darken_cmap(plt.cm.RdBu_r, 0.9)
        else:
            cmap_use = color_map
    if order_points:
        o = np.argsort(color_dat)
    else:
        o = np.arange(len(color_dat)) 
    pp = ax.scatter(x[o], y[o], s=marker_size, c = color_dat[o], 
        cmap=cmap_use, vmin=vmin, vmax=vmax)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Doublet score')
    ax.set_xlabel(embedding_name + ' 1')
    ax.set_ylabel(embedding_name + ' 2')
    fig.colorbar(pp, ax=ax)

    ax = axs[0]
    called_doubs = self.predicted_doublets_
    ax.scatter(x[o], y[o], s=marker_size,  c=called_doubs[o], cmap=custom_cmap([[.7,.7,.7], [0,0,0]]))
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Predicted doublets')
    #singlet_marker = Line2D([], [], color=[.7,.7,.7], marker='o', markersize=5, label='Singlet', linewidth=0)
    #doublet_marker = Line2D([], [], color=[.0,.0,.0], marker='o', markersize=5, label='Doublet', linewidth=0)
    #ax.legend(handles = [singlet_marker, doublet_marker])
    ax.set_xlabel(embedding_name + ' 1')
    ax.set_ylabel(embedding_name + ' 2')


# Load data
counts_path = "Filtered_peaks_bc_matrixes/{}/filtered_peak_bc_matrix/matrix.mtx".format(gem_id)
barcodes_path = "Filtered_peaks_bc_matrixes/{}/filtered_peak_bc_matrix/barcodes.tsv".format(gem_id) 
counts_matrix = scipy.io.mmread(counts_path).T.tocsc()
barcodes_df = pd.read_csv(barcodes_path, header = None)

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    

if gem_id == "m7maqeew_gt2mrtf3" or gem_id == "xazsjcvt_kzi9e2rf":
    expected_doublet_rate = 0.056

else:
    expected_doublet_rate = 0.04


scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(log_transform=True, min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=70, 
                                                          n_prin_comps=50)
doublet_scores = np.round(doublet_scores, decimals = 3)


# Create data frame
scrublet_doubl_dict = {"barcodes":barcodes_df[0].values, "scrublet_doublet_scores":doublet_scores, "scrublet_predicted_doublet":predicted_doublets}
scrublet_doubl_df = pd.DataFrame(scrublet_doubl_dict)

# Save
if not os.path.exists("./tmp"):    
    os.mkdir("./tmp")
if not os.path.exists("./tmp/histograms/"):
    os.mkdir("./tmp/histograms/")
if not os.path.exists("./tmp/umaps/"):
    os.mkdir("./tmp/umaps/")

# Plot
scrub.plot_histogram()
plt.savefig("tmp/histograms/scrublet_doublet_prediction_histograms-{}.png".format(gem_id), dpi = 100)

scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
plot_embedding(scrub, 'UMAP', order_points=True)
plt.savefig("tmp/umaps/scrublet_doublet_prediction_umaps-{}.png".format(gem_id), dpi = 100)

scrublet_doubl_df.to_csv("tables/scrublet_doublet_prediction-{}.csv".format(gem_id), index = False)

ending_time = datetime.datetime.now()

print("Ending job at " + ending_time.strftime("%Y-%m-%d %H:%M:%S"))



