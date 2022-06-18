#####################################################################################
########################## Doublet Detection with Scrublet ##########################
#####################################################################################

# Description:
# This script reads the filtered_feature_bc_matrix from a specific GEM-id and runs scrublet.
# It saves a dataframe with the doublet score and prediction for each cell, the UMAPs with the
# scores projected and the histograms of their distribution.


# Import packages
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


# Initialize variables
subproject = sys.argv[1]
gem_id = sys.argv[2]


# Load data
counts_path = "../../1-cellranger_mapping/projects/{}/jobs/{}/{}/outs/filtered_feature_bc_matrix/matrix.mtx.gz".format(subproject, gem_id, gem_id)
barcodes_path = "../../1-cellranger_mapping/projects/{}/jobs/{}/{}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz".format(subproject, gem_id, gem_id) 
counts_matrix = scipy.io.mmread(counts_path).T.tocsc()
barcodes_df = pd.read_csv(barcodes_path, header = None)
metadata = pd.read_csv("../../1-cellranger_mapping/data/tonsil_atlas_metadata.csv")


# Run scrublet
lib_type = metadata.loc[metadata["gem_id"] == gem_id, :]["type"].values[0]
if lib_type == "not_hashed":
	expected_doublet_rate = 0.04
elif lib_type == "hashed_cdna" or lib_type == "hashed_hto":
	expected_doublet_rate = 0.16
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate = expected_doublet_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 75, n_prin_comps = 30)
doublet_scores = np.round(doublet_scores, decimals = 3)


# Get 2D embedding to visualize results
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))


# Create data frame
scrublet_doubl_dict = {"barcodes":barcodes_df[0].values, "scrublet_doublet_scores":doublet_scores, "scrublet_predicted_doublet":predicted_doublets}
scrublet_doubl_df = pd.DataFrame(scrublet_doubl_dict)


# Plot
hists = scrub.plot_histogram()
umaps = scrub.plot_embedding('UMAP', order_points=True)


# Save
if not os.path.exists("./tmp"):    
    os.mkdir("./tmp")
if not os.path.exists("./tmp/histograms/"):
    os.mkdir("./tmp/histograms/")
if not os.path.exists("./tmp/umaps/"):
    os.mkdir("./tmp/umaps/")
if not os.path.exists("../../results/tables/scrublet"):
    os.mkdir("../../results/tables/scrublet/")
scrublet_doubl_df.to_csv("../../results/tables/scrublet/scrublet_doublet_prediction-{}-{}.csv".format(subproject, gem_id), index = False)
hists[0].savefig("./tmp/histograms/scrublet_doublet_prediction_histograms-{}-{}.png".format(subproject, gem_id), dpi = 100)
umaps[0].savefig("./tmp/umaps/scrublet_doublet_prediction_umaps-{}-{}.png".format(subproject, gem_id), dpi = 100)


ending_time = datetime.datetime.now()
print("Ending job at " + ending_time.strftime("%Y-%m-%d %H:%M:%S"))
