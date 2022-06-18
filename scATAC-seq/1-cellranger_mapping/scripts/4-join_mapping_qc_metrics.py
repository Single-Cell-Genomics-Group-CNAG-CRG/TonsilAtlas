# This script joins all the summary.csv of all libraries in a single file


# Load modules
import pandas as pd
import numpy as np
import os


# Load data
metadata = pd.read_csv("data/tonsil_atlas_metadata_atac.csv")
iterable = zip(metadata["subproject"], metadata["gem_id"])
summary_path = "projects/{}/jobs/{}/{}/outs/summary.csv"
summary_files_paths = [summary_path.format(subproject, gem_id, gem_id) for subproject, gem_id in iterable]
summary_files_dfs = [pd.read_csv(path) for path in summary_files_paths]


# Merge
summary_df = summary_files_dfs[0]
for i in range(1, len(summary_files_dfs)):
    summary_df = summary_df.append(summary_files_dfs[i])


# Save
gem_id_col = list(metadata["gem_id"])
summary_df.insert(0, column = "gem_id", value = gem_id_col)
summary_df.to_csv("../results/tables/cellranger_mapping/cellranger_mapping_metrics_atac.csv", header = True, index = None)
