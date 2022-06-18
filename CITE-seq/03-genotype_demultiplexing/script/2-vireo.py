"""
Author: "Sonal Rashmi"
Date: 2021-05-18
Description: This script takes the CellSNP output directory and metadata as input, based on the SNP profile demultiplex the cells for each sample in metadata. It further classifies cells as doublet and unassigned based on the genotype profile.  
"""


#!/usr/bin/python

import sys
import vireoSNP
import numpy as np
from scipy import sparse
from scipy.io import mmread
import pandas as pd 
import os
import argparse
import matplotlib.pyplot as plt
from vireoSNP.plot.base_plot import heat_matrix
		

class VireoDonorDeconvolution(object):
	def __init__(self, cellsnp_path, gem_id, vireo_output_dir, donor_n = None):
		self.ad_matrix=os.path.join(cellsnp_path,gem_id,"cellSNP.tag.AD.mtx")
		self.dp_matrix=os.path.join(cellsnp_path,gem_id,"cellSNP.tag.DP.mtx")
		self.vcf=os.path.join(cellsnp_path,gem_id,"cellSNP.vcf.gz")
		self.outdir = vireo_output_dir
		self.gem_id = gem_id
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)
		
		if os.path.exists(self.ad_matrix) and os.path.exists(self.dp_matrix):
			try:		
				self.cell_vcf = vireoSNP.load_VCF(self.vcf, biallelic_only=True)
				self.cell_dat = vireoSNP.vcf.read_sparse_GeneINFO(cell_vcf['GenoINFO'], keys=['AD', 'DP'])

				self.AD = cell_dat['AD']
				self.DP = cell_dat['DP']
			except:	
				self.AD = mmread(self.ad_matrix).tocsc()
				self.DP = mmread(self.dp_matrix).tocsc()
			if donor_n == None:
				self.donor_n = 4
			self.res = vireoSNP.vireo_wrap(self.AD, self.DP, n_donor=self.donor_n, learn_GT=True,
                          n_extra_donor=0, ASE_mode=False, fix_beta_sum=False,
                          n_init=50, check_doublet=True, random_seed=1)
			cmd="vireo -c "+os.path.join(cellsnp_path,gem_id)+" -N 4 -o "+os.path.join(vireo_output_dir,gem_id)+" --randSeed 2"                          	
			os.system(cmd)
		else:
			print("CellSNP files not found")
			sys.exit(0)	

	def get_donor_number(self):
			
		n_donor_list = np.arange(2, 6)
		ELBO_list_all = []
		for _n_don in n_donor_list:
			res = vireoSNP.vireo_wrap(self.AD, self.DP, n_donor=_n_don, learn_GT=True,
									  n_extra_donor=0, ASE_mode=False, fix_beta_sum=False,
									  n_init=50, check_doublet=True, random_seed=1)
			ELBO_list_all.append(res['LB_list'])
			
			
		fig = plt.figure(figsize=(5, 4), dpi=100)
		plt.plot(n_donor_list - 1, np.max(ELBO_list_all, axis=1))
		plt.boxplot(ELBO_list_all)
		plt.xticks(n_donor_list - 1, n_donor_list)
		plt.ylabel("ELBO")
		plt.xlabel("n_clones")
		plt.savefig(os.path.join(self.outdir,(self.gem_id+".pdf")))
		
	
	def donor_assignment_probability(self):
		
		fig = plt.figure(figsize=(4, 5), dpi=100)
		# assign_prob_comb = self.res['ID_prob']
		assign_prob_comb = np.append(self.res['ID_prob'], 
									 np.sum(self.res['doublet_prob'], axis=1, 
											keepdims=True), axis=1)
		im = heat_matrix(assign_prob_comb, 
						 cmap="Oranges", alpha=0.8,
						 display_value=False, row_sort=True)
		plt.colorbar(im, fraction=0.046, pad=0.04)
		plt.title("cell assignment prob")
		plt.xlabel("Donors")
		plt.ylabel("%d cells" %(self.res['ID_prob'].shape[0]))
		plt.yticks([])
		plt.xticks([0, 1, 2, 3, 4], [1, 2, 3, 4, "doublet"])

		# plt.tight_layout()
		plt.savefig(os.path.join(self.outdir,(self.gem_id+"_cell_assignment.pdf")))
	
	def get_number_of_doublet(self):
		doublet_threshold = 0.9
		is_doublet = np.sum(self.res['doublet_prob'], axis=1) > doublet_threshold
		print("%d cells are called doublet" %(sum(is_doublet)))
		return sum(is_doublet)

	def get_number_of_cells_unassigned(self):
		prob_threshold = 0.9
		is_doublet = np.sum(self.res['doublet_prob'], axis=1) > prob_threshold   
		is_unassigned = (np.max(self.res['ID_prob'], axis=1) < prob_threshold) & (~is_doublet)
		print("%d cells are unassigned to singlet or doublets" %(sum(is_unassigned)))
		return sum(is_unassigned)
		
	def get_allelic_ratio_per_variant_per_donor(self):
		## If ASE_mode is False
		AF_SNPs = np.tensordot(self.res['GT_prob'], self.res['theta_mean'][0, :], axes=[2, 0])
		
		fig = plt.figure(figsize=(4, 5), dpi=100)
		im = heat_matrix(AF_SNPs, cmap="GnBu", alpha=0.8,
						 display_value=False, row_sort=True)
		plt.colorbar(im, fraction=0.046, pad=0.04)
		plt.title("Mean allelic ratio")
		plt.xlabel("Donors")
		plt.ylabel("%d SNPs" %(AF_SNPs.shape[0]))
		plt.yticks([])
		plt.xticks([0, 1, 2, 3])

		plt.tight_layout()
		plt.savefig(os.path.join(self.outdir,(self.gem_id+"_heatmap_profile.pdf")))
	
	
def main():
	parser = argparse.ArgumentParser(description = "options to run vireo using cellsnp")
	parser.add_argument("--cellsnp_output_path",
						dest = "cellsnp_output_path",
						action = "store",
						default = None,
						help = "Successful run cellsnp output")
	parser.add_argument("--metadata_path",
						dest = "metadata_path",
						action = "store",
						default = None,
						help = "Metadata file used to run cellranger")					
	parser.add_argument("--vireo_output_dir",
						dest = "vireo_output_dir",
						action = "store",
						default = None,
						help = "Vireo Output Dir")
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(0)
	options = parser.parse_args()
	metadata_path = options.metadata_path
	cellsnp_output_path = options.cellsnp_output_path
	vireo_output_dir = options.vireo_output_dir	
	if not os.path.exists(metadata_path):
		print("metadata file not found. Program terminated!")
		sys.exit(0)
	else:
		metadata = pd.read_csv(metadata_path)
		metadata_multiplexed = metadata[metadata['donor_id'] == "multiplexed"]
		iterable = list(sorted(set(zip(metadata_multiplexed["subproject"], metadata_multiplexed["gem_id"]))))
		d=[]   
		for subproject, gem_id in iterable:
			VireoDonorDeconvolutionObj=VireoDonorDeconvolution(cellsnp_output_path, gem_id, vireo_output_dir)
			VireoDonorDeconvolutionObj.get_donor_number()
			VireoDonorDeconvolutionObj.donor_assignment_probability()
			doublet = VireoDonorDeconvolutionObj.get_number_of_doublet()
			unassigned = VireoDonorDeconvolutionObj.get_number_of_cells_unassigned()
			#df=pd.DataFrame(VireoDonorDeconvolutionObj.res)
			#df.to_csv(os.path.join(vireo_output_dir,(gem_id+"_df.csv")))
			d.append({'Gem_id': gem_id, 'Doublet': doublet,'Unassigned':  unassigned})
		summary_df = pd.DataFrame(d)
		summary_df.to_csv(os.path.join(vireo_output_dir,"Vireo_Summary.csv"))    			
      
if __name__ =="__main__":
		main()            
   