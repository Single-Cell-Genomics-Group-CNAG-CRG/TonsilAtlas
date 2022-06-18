# This script joins all the summary.csv of all libraries in a single file

"""
Author: "Sonal Rashmi"
Date: 2021-05-12
Description: This script merges metrics of multi cellranger 6.0.1 for all libraries into a single file
Input Files : Metadata file and cellranger working directory
Output Files: cellranger_mapping_metrics_(counts/vdj_t/vdj_b/adt).csv in results (creates results folder if not present)
"""

# Load modules
import pandas as pd
import numpy as np
import os
import argparse
import sys

def main():
	parser = argparse.ArgumentParser(description='The script merges metrics of multi cellranger 6.0.1 for all gem_ids into respective lib_type files')
	parser.add_argument('-metadata_file',const=None, help="metadata_file containing gem_id and subproject", dest='inputmetadatafile')
	parser.add_argument('-working_directory',const=None, help="Cellranger output directory with structure as projects/subproject/jobs/gem_id", dest='workdir')
	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(0)
	args=parser.parse_args()
	input_metadata_file = args.inputmetadatafile
	work_dir = args.workdir
	if os.path.exists(input_metadata_file) and os.path.exists(work_dir):
		# Load data
		metadata = pd.read_csv(input_metadata_file)
		iterable = list(sorted(set(zip(metadata["subproject"], metadata["gem_id"]))))

		summary_gex_dfs = pd.DataFrame()
		summary_vdj_t_dfs = pd.DataFrame()
		summary_vdj_b_dfs = pd.DataFrame()
		summary_adt_dfs = pd.DataFrame()
		for subproject, gem_id in iterable:
			summary_path = os.path.join(work_dir,"projects/{}/jobs/{}/{}/outs/per_sample_outs/{}/metrics_summary.csv".format(subproject, gem_id, gem_id, gem_id))
			if os.path.exists(summary_path):
				dat = pd.read_csv(summary_path)
				dat.columns = dat.columns.str.replace(' ','_')
				dat['Metric_Name']=dat['Library_or_Sample']+'_'+dat['Metric_Name']
				for i in dat.columns:
					dat[i] = dat[i].str.replace(' ','_')
				lib_type = dat['Library_Type'].to_list()	
				if 'VDJ_T' in lib_type:
					vdj_t_df = dat[(dat['Library_or_Sample'] == 'Library') & (dat['Library_Type'] == 'VDJ_T') & (dat['Grouped_By'] == 'Physical_library_ID') | (dat['Library_or_Sample'] == 'Sample') & (dat['Library_Type'] == 'VDJ_T') ][['Metric_Name','Metric_Value']]
					vdj_t_df = vdj_t_df.set_index('Metric_Name')
					vdj_t_df.index.names = [None]
					vdj_t_df.columns = [gem_id]
					vdj_t_df_t = vdj_t_df.transpose()
					vdj_t_df_t.reset_index(inplace=True)
					vdj_t_df_t = vdj_t_df_t.rename(columns = {'index':'gem_id'})
					vdj_t_df_t.insert(0, column = "subproject", value = subproject)
					summary_vdj_t_dfs = summary_vdj_t_dfs.append(vdj_t_df_t)
				if 'VDJ_B' in lib_type:
					vdj_b_df = dat[(dat['Library_or_Sample'] == 'Library') & (dat['Library_Type'] == 'VDJ_B') & (dat['Grouped_By'] == 'Physical_library_ID') | (dat['Library_or_Sample'] == 'Sample') & (dat['Library_Type'] == 'VDJ_B') ][['Metric_Name','Metric_Value']]
					vdj_b_df = vdj_b_df.set_index('Metric_Name')
					vdj_b_df.index.names = [None]
					vdj_b_df.columns = [gem_id]
					vdj_b_df_t = vdj_b_df.transpose()
					vdj_b_df_t.reset_index(inplace=True)
					vdj_b_df_t = vdj_b_df_t.rename(columns = {'index':'gem_id'})
					vdj_b_df_t.insert(0, column = "subproject", value = subproject)
					summary_vdj_b_dfs = summary_vdj_b_dfs.append(vdj_b_df_t)
				if 'Gene_Expression' in lib_type:
					gex_df = dat[(dat['Library_or_Sample'] == 'Library') & (dat['Library_Type'] == 'Gene_Expression') & (dat['Grouped_By'] == 'Physical_library_ID') | (dat['Library_or_Sample'] == 'Sample') & (dat['Library_Type'] == 'Gene_Expression') ][['Metric_Name','Metric_Value']]	
					gex_df = gex_df.set_index('Metric_Name')
					gex_df.index.names = [None]
					gex_df.columns = [gem_id]
					gex_df_t = gex_df.transpose()
					gex_df_t.reset_index(inplace=True)
					gex_df_t = gex_df_t.rename(columns = {'index':'gem_id'})
					gex_df_t.insert(0, column = "subproject", value = subproject)
					summary_gex_dfs = summary_gex_dfs.append(gex_df_t)
				if 'Antibody_Capture' in lib_type:
					adt_df = dat[(dat['Library_or_Sample'] == 'Library') & (dat['Library_Type'] == 'Antibody_Capture') & (dat['Grouped_By'] == 'Physical_library_ID') | (dat['Library_or_Sample'] == 'Sample') & (dat['Library_Type'] == 'Antibody_Capture') ][['Metric_Name','Metric_Value']]
					adt_df = adt_df.set_index('Metric_Name')
					adt_df.index.names = [None]
					adt_df.columns = [gem_id]
					adt_df_t = adt_df.transpose()
					adt_df_t.reset_index(inplace=True)
					adt_df_t = adt_df_t.rename(columns = {'index':'gem_id'})
					adt_df_t.insert(0, column = "subproject", value = subproject)
					summary_adt_dfs = summary_adt_dfs.append(adt_df_t)
			else:
				print("Project Job folders not found for the gemid "+ gem_id) 
		if not os.path.exists(os.path.join(work_dir,'results/')):
			os.makedirs(os.path.join(work_dir,'results/'))
		if not summary_gex_dfs.empty:
			summary_gex_dfs.to_csv(os.path.join(work_dir,'results',"cellranger_mapping_metrics_gex.csv"), header = True, index = None)
		if not summary_vdj_b_dfs.empty:	
			summary_vdj_b_dfs.to_csv(os.path.join(work_dir,'results',"cellranger_mapping_metrics_vdj_b.csv"), header = True, index = None)
		if not summary_vdj_t_dfs.empty:
			summary_vdj_t_dfs.to_csv(os.path.join(work_dir,'results',"cellranger_mapping_metrics_vdj_t.csv"), header = True, index = None)
		if not summary_adt_dfs.empty:
			summary_adt_dfs.to_csv(os.path.join(work_dir,'results',"cellranger_mapping_metrics_adt.csv"), header = True, index = None)

if __name__ =="__main__":
		main()			
			