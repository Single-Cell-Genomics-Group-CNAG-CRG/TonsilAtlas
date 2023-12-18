# This script clones the TonsilAtlas github repo and generates a csv file
# with all the fastq files that we need to upload to HCA


# Load relevant packages
import os
import numpy as np
import pandas as pd
from git import Repo  # pip install gitpython


# Define variables
ordered_cols = ["dataset", "subproject", "gem_id", "library_id",
                "library_name", "type_x", "pair_id", "read", "fastq_path",
                "donor_id"]


# Clone TonsilAtlas repo
git_url = "https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/"
repo_dir = os.getcwd()
Repo.clone_from(git_url, repo_dir)


# Read all fastq paths: RNA
path_rna_metadata = f"{repo_dir}/scRNA-seq/1-cellranger_mapping/data/tonsil_atlas_metadata.csv"
rna_metadata = pd.read_csv(path_rna_metadata)
subprojects = np.unique(rna_metadata["subproject"])
fastqs_rna_list = []
for subproject in subprojects:
    path_fastqs_csv = f"{repo_dir}/scRNA-seq/1-cellranger_mapping/projects/{subproject}/fastq_paths.csv"
    fastq_csv = pd.read_csv(path_fastqs_csv)
    fastqs_rna_list.append(fastq_csv)

fastqs_rna = pd.concat(fastqs_rna_list) 
fastqs_rna = pd.merge(fastqs_rna, rna_metadata, on = "library_id", how = "left")
fastqs_rna["dataset"] = "scRNA-seq"
fastqs_rna = fastqs_rna[ordered_cols]
fastqs_rna.rename(columns = {"type_x": "type"}, inplace = True)


# Read all fastq paths: scATAC-seq
path_atac_metadata = f"{repo_dir}/scATAC-seq/1-cellranger_mapping/data/tonsil_atlas_metadata_atac.csv"
atac_metadata = pd.read_csv(path_atac_metadata)
subprojects = np.unique(atac_metadata["subproject"])
fastqs_atac_list = []
for subproject in subprojects:
    path_fastqs_csv = f"{repo_dir}/scATAC-seq/1-cellranger_mapping/projects/{subproject}/fastq_paths.csv"
    fastq_csv = pd.read_csv(path_fastqs_csv)
    fastqs_atac_list.append(fastq_csv)

fastqs_atac = pd.concat(fastqs_atac_list) 
fastqs_atac = pd.merge(fastqs_atac, atac_metadata, on = "library_id", how = "left")
fastqs_atac["dataset"] = "scATAC-seq"
fastqs_atac["type_x"] = "ATAC"
fastqs_atac = fastqs_atac[ordered_cols]
fastqs_atac.rename(columns = {"type_x": "type"}, inplace = True)


# Multiome
path_multiome_metadata = f"{repo_dir}/multiome/1-cellranger_mapping/data/tonsil_atlas_metadata_multiome.csv"
multiome_metadata = pd.read_csv(path_multiome_metadata)
subprojects = np.unique(multiome_metadata["subproject"])
multiome_dict = {
    "BCLLATLAS_133":f"{repo_dir}/multiome/1-cellranger_mapping/projects/experiment_3/fastq_paths_rna.csv",
    "BCLLATLAS_134":f"{repo_dir}/multiome/1-cellranger_mapping/projects/experiment_3/fastq_paths_atac.csv",
    "BCLLATLAS_42":f"{repo_dir}/multiome/1-cellranger_mapping/projects/experiment_1/fastq_paths_rna.csv",
    "BCLLATLAS_43":f"{repo_dir}/multiome/1-cellranger_mapping/projects/experiment_1/fastq_paths_atac.csv",
    "BCLLATLAS_47":"multiome/1-cellranger_mapping/projects/experiment_2/fastq_paths_rna.csv",
    "BCLLATLAS_48":"multiome/1-cellranger_mapping/projects/experiment_2/fastq_paths_atac.csv",
}
fastqs_multiome_list = []
for subproject in subprojects:
    fastq_csv = pd.read_csv(multiome_dict[subproject])
    fastqs_multiome_list.append(fastq_csv)

fastqs_multiome = pd.concat(fastqs_multiome_list) 
fastqs_multiome = pd.merge(fastqs_multiome, multiome_metadata, on = "library_id", how = "left")
fastqs_multiome["dataset"] = "multiome"
fastqs_multiome.rename(columns = {"type": "type_x"}, inplace = True)
fastqs_multiome = fastqs_multiome[ordered_cols]
fastqs_multiome.rename(columns = {"type_x": "type"}, inplace = True)


# CITE-seq
path_cite_metadata = f"{repo_dir}/CITE-seq/01-cellranger_mapping/data/tonsil_atlas_metadata.csv"
cite_metadata = pd.read_csv(path_cite_metadata)
subprojects = np.unique(cite_metadata["subproject"])
fastqs_cite_list = []
for subproject in subprojects:
    path_fastqs_csv = f"{repo_dir}/CITE-seq/01-cellranger_mapping/projects/{subproject}/fastq_paths.csv"
    fastq_csv = pd.read_csv(path_fastqs_csv)
    fastqs_cite_list.append(fastq_csv)

fastqs_cite = pd.concat(fastqs_cite_list) 
fastqs_cite = pd.merge(fastqs_cite, cite_metadata, on = "library_id", how = "left")
fastqs_cite["dataset"] = "CITE-seq"
fastqs_cite = fastqs_cite[ordered_cols]
fastqs_cite.rename(columns = {"type_x": "type"}, inplace = True)


# MCL
path_mcl_metadata = f"{repo_dir}/MCL/1-cellranger_mapping/data/sequencing_metadata.csv"
mcl_metadata = pd.read_csv(path_mcl_metadata)
subprojects = np.unique(mcl_metadata["subproject"])
mcl_dict = {
    "BCLLATLAS_64":f"{repo_dir}/MCL/1-cellranger_mapping/projects/experiment_1/fastq_paths_rna.csv",
    "BCLLATLAS_65":f"{repo_dir}/MCL/1-cellranger_mapping/projects/experiment_1/fastq_paths_atac.csv",
}
fastqs_mcl_list = []
for subproject in subprojects:
    fastq_csv = pd.read_csv(mcl_dict[subproject])
    fastqs_mcl_list.append(fastq_csv)

fastqs_mcl = pd.concat(fastqs_mcl_list) 
fastqs_mcl = pd.merge(fastqs_mcl, mcl_metadata, on = "library_id", how = "left")
fastqs_mcl["dataset"] = "MCL"
fastqs_mcl.rename(columns = {"type": "type_x"}, inplace = True)
fastqs_mcl = fastqs_mcl[ordered_cols]
fastqs_mcl.rename(columns = {"type_x": "type"}, inplace = True)


# Spatial
path_spatial_metadata = f"{repo_dir}/spatial_transcriptomics/01-spaceranger/data/sample_id.txt"
spatial_metadata = pd.read_csv(path_spatial_metadata)
path_fastqs_spatial = f"{repo_dir}/spatial_transcriptomics/01-spaceranger/data/fastq_paths.csv"
fastqs_spatial = pd.read_csv(path_fastqs_spatial)
fastqs_spatial.rename(columns = {"type": "type_x", "gemid":"gem_id"}, inplace = True)
fastqs_spatial["dataset"] = "spatial_transcriptomics"
spatial_metadata = spatial_metadata[["gem_id", "donor_id"]]
fastqs_spatial = pd.merge(fastqs_spatial, spatial_metadata, on = "gem_id", how = "left")
fastqs_spatial["library_name"] = [f"{x}_spatial" for x in fastqs_spatial.donor_id]
def extract_identifier(path):
    return path.rsplit('_', 1)[0]

fastqs_spatial['identifier'] = fastqs_spatial['fastq_path'].apply(extract_identifier)
unique_identifiers = fastqs_spatial['identifier'].unique()
pair_id_map = {id: f"P{i+1}" for i, id in enumerate(unique_identifiers)}
fastqs_spatial['pair_id'] = fastqs_spatial['identifier'].map(pair_id_map)
fastqs_spatial.drop('identifier', axis=1, inplace=True)
fastqs_spatial = fastqs_spatial[ordered_cols]
fastqs_spatial.rename(columns = {"type_x": "type"}, inplace = True)


# Concat all
fastq_all = pd.concat([fastqs_rna, fastqs_atac, fastqs_multiome, fastqs_cite, fastqs_spatial, fastqs_mcl])


# Specify flowcell, lane and index
flowcell = []
lane = []
index = []
file_paths = fastq_all["fastq_path"]
for path in file_paths:
    parts = path.split('/')
    file_name = parts[-1]
    file_parts = file_name.split('_')
    flowcell.append(file_parts[0])
    lane.append(file_parts[1].split('-')[0])
    index.append(file_parts[2])

fastq_all["flowcell"] = flowcell
fastq_all["lane"] = lane
fastq_all["index"] = index


# Include name to upload to ArrayExpress
scRNAseq.BCLL-2-T.BCLLATLAS_05.jb6vuao4_g4vi9ur0.130339.NotHashed.P1.R1.fastq.gz
fastq_name_array_express = []
for i in range(fastq_all.shape[0]):
    dataset = fastq_all["dataset"].values[i]
    donor_id = fastq_all["donor_id"].values[i]
    subproject = fastq_all["subproject"].values[i]
    gem_id = fastq_all["gem_id"].values[i]
    library_id = fastq_all["library_id"].values[i]
    lib_type = fastq_all["type"].values[i]
    pair_id = fastq_all["pair_id"].values[i]
    read = fastq_all["read"].values[i]
    tmp = f"{dataset}.{donor_id}.{subproject}.{gem_id}.{library_id}.{lib_type}.{pair_id}.{read}.fastq.gz"
    fastq_name_array_express.append(tmp)

fastq_all["fastq_name_array_express"] = fastq_name_array_express


# Rename read 1 and 2 from spatial transcriptomics
fastq_all.read[fastq_all.read == 1] = "R1"
fastq_all.read[fastq_all.read == 2] = "R2"

# Save
fastq_all.to_csv(f"{repo_dir}/fastq_paths_all.csv", index = False)
