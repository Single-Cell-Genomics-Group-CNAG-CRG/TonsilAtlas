#####################################################################################
########################## Doublet Detection with Scrublet ##########################
#####################################################################################

# Description:
# This script reads the sparse matrices for RNA and ATAC of the multiome samples and runs scrublet.
# It saves a dataframe with the doublet score and prediction for each cell, the UMAPs with the
# scores projected and the histograms of their distribution.


# Import packages
print("Importing packages...")
import numpy as np
import pandas as pd
import scipy.io
import sys
import os
import scrublet as scr
import pynndescent
import datetime


starting_time = datetime.datetime.now()
print("Starting job at " + starting_time.strftime("%Y-%m-%d %H:%M:%S"))


# Load data
print("Loading data...")
all_data = {}
all_keys = os.listdir("tmp")
for key in all_keys:
    print(key)
    all_data[key] = []
    path_to_matrix = "tmp/{}/matrix.mtx".format(key)
    path_to_barcodes = "tmp/{}/barcodes.tsv".format(key)
    path_to_features = "tmp/{}/features.tsv".format(key)
    matrix = scipy.io.mmread(path_to_matrix).T.tocsc()
    barcodes = pd.read_csv(path_to_barcodes, header = None)
    features = pd.read_csv(path_to_features, header = None)
    all_data[key].extend([matrix, barcodes, features])



# Create directories
print("Creating directories...")
if not os.path.exists("../../results/tables/scrublet"):
    os.mkdir("../../results/tables/scrublet/")
if not os.path.exists("./tmp/histograms/"):
    os.mkdir("./tmp/histograms/")
if not os.path.exists("./tmp/umaps/"):
    os.mkdir("./tmp/umaps/")


# Run
print("Running scrublet...")
expected_doublet_rate = 0.056
for key in all_data.keys():
    print(key)
    # Scrublet
    scrub = scr.Scrublet(all_data[key][0], expected_doublet_rate = expected_doublet_rate)
    if "atac" in key:
        doublet_scores, predicted_doublets = scrub.scrub_doublets(log_transform = True, min_counts = 2, min_cells = 3, min_gene_variability_pctl = 70, n_prin_comps = 50)
    elif "rna" in key:
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 75, n_prin_comps = 30)
    doublet_scores = np.round(doublet_scores, decimals = 3)
    
    
    # Create DataFrame
    scrublet_doubl_dict = {"barcodes": all_data[key][1][0].values, "scrublet_doublet_scores": doublet_scores, "scrublet_predicted_doublet": predicted_doublets}
    scrublet_doubl_df = pd.DataFrame(scrublet_doubl_dict)
    all_data[key].append(scrublet_doubl_df)
    
    
    # Plot
    scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist = 0.3))
    hists = scrub.plot_histogram()
    umaps = scrub.plot_embedding("UMAP", order_points = True)
    

    # Save
    scrublet_doubl_df.to_csv("../../results/tables/scrublet/scrublet_doublet_prediction_{}.csv".format(key), index = False)
    hists[0].savefig("./tmp/histograms/scrublet_doublet_prediction_histograms_{}.png".format(key), dpi = 100)
    umaps[0].savefig("./tmp/umaps/scrublet_doublet_prediction_umaps_{}".format(key), dpi = 100)


ending_time = datetime.datetime.now()
print("Ending job at " + ending_time.strftime("%Y-%m-%d %H:%M:%S"))
