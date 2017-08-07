#########################################################################################
#																						#		
#				Network Building stats													#
#																						#
#										Jim Brunner 									#
#																						#
#																						#
#																						#
#########################################################################################

## modules needed
from pylab import *
import pandas as pd
import sys
import matplotlib.pyplot as plt
from scipy import misc, sparse, cluster
import os
from co_occ_funs import *
from random import sample
import time
from joblib import Parallel, delayed, cpu_count

if __name__ == '__main__':

	#import networkx as nx

	#####Variables
	##User input, should be entered as an argument
	csv_name = sys.argv[1]#the data
	level = 'species'


	######Import the abundance matrix
	#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 
	abundance_array_full = pd.read_csv(csv_name, sep = ' ')
	#Let's remove taxa that aren't present at all.
	for ind in abundance_array_full.index:
		if max(abundance_array_full.loc[ind][2:]) == 0:
			abundance_array_full = abundance_array_full.drop(ind)

	abundance_array_full.index = range(len(abundance_array_full))		


	the_levels = abundance_array_full['LEVEL'].values
	the_taxa = abundance_array_full['TAXA'].values


		
	num_tris = 20
	


	abundance_array = abundance_array_full
	#start by getting the indices of the members of the level
	in_level = where(abundance_array['LEVEL'] == level)[0]
	#then the part of the dataframe containing just that level
	lvl_abundance_array = abundance_array.iloc[in_level]
	#create np array of abundance values (dropping labels) within the level
	ab_np_array1 = lvl_abundance_array.values[:,2:]
	not_seen = where(sum(ab_np_array1,1) == 0)[0]
	n_seen_ind = lvl_abundance_array.index[not_seen]
	in_level = delete(in_level,not_seen)
	#prune off unseen organisms
	lvl_abundance_array = lvl_abundance_array.drop(n_seen_ind)

	num_samples = len(lvl_abundance_array.columns) - 2
	num_tax = len(lvl_abundance_array)
	print(num_tax)

	num_edges_p = zeros([num_tris,num_samples])
	mean_deg_p = zeros([num_tris,num_samples])
	p_md_p = zeros([num_tris,num_samples])
	p_vd_p = zeros([num_tris,num_samples])
	p_max_ev_p = zeros([num_tris,num_samples])
	p_min_ev_p = zeros([num_tris,num_samples])
	p_num_edges_p = zeros([num_tris,num_samples])

	num_edges_pt = zeros([num_tris,num_samples])
	mean_deg_pt = zeros([num_tris,num_samples])
	p_md_pt = zeros([num_tris,num_samples])
	p_vd_pt = zeros([num_tris,num_samples])
	p_max_ev_pt = zeros([num_tris,num_samples])
	p_min_ev_pt = zeros([num_tris,num_samples])
	p_num_edges_pt = zeros([num_tris,num_samples])

	# 		num_edges_b = zeros(num_samples)
	# 		mean_deg_b = zeros(num_samples)
	# 		p_md_b = zeros(num_samples)
	# 		p_vd_b = zeros(num_samples)
	# 		p_max_ev_b = zeros(num_samples)
	# 		p_min_ev_b = zeros(num_samples)	
	#		p_num_edges_b = zeros(num_samples)
	elapseds = zeros([num_tris,num_samples])

	for n in range(num_tris):
		for ns in range(2,num_samples):
			print('\r',ns, '/',num_samples - 1, end = '     ')
			st = time.perf_counter()
			chs = array(sample(range(num_samples), ns)) + 2
			this_array = lvl_abundance_array.loc[:,lvl_abundance_array.columns[append([0,1],chs)]]
			#print(this_array.columns)
			#build a network with a subset of the data
			pears = build_network(this_array, 'pearson', thr = False, list_too = False)
			pears_thr = build_network(this_array, 'pearson', thr = True, list_too = False)
			#binned = build_network(this_array, 'bins', list_too = False)

			pears_data = pears.values
			pears_thr_data = pears_thr.values
			#binned_data = binned.values

			num_edges_p[n,ns] = len(pears_data.nonzero()[0])/2
			num_edges_pt[n,ns] = len(pears_thr_data.nonzero()[0])/2
			#num_edges_b[n,ns] = len(binned_data.nonzero()[0])/2

			mean_deg_p[n,ns] = mean(pears_data.sum(axis = 1))
			mean_deg_pt[n,ns] = mean(pears_thr_data.sum(axis = 1))
			#mean_deg_b[n,ns] = mean(binned_data.sum(axis = 1))

			num_mc_sims = 1000
			this_np_array = this_array.values[:,2:]
			[p_md_p[n,ns], p_vd_p[n,ns], p_max_ev_p[n,ns], p_min_ev_p[n,ns],p_num_edges_p[n,ns]] = mc_network_stats(this_np_array, pears_data, sims = num_mc_sims)
			[p_md_pt[n,ns], p_vd_pt[n,ns], p_max_ev_pt[n,ns], p_min_ev_pt[n,ns], p_num_edges_pt[n,ns]] = mc_network_stats(this_np_array, pears_thr_data, thr = True, sims = num_mc_sims)
			#[p_md_b[n,ns], p_vd_b[n,ns], p_max_ev_b[n,ns], p_min_ev_b[n,ns],p_num_edges_b[n,ns]] = mc_network_stats(this_np_array, binned_data, bins = True, sims = num_mc_sims)
	
			elapseds[n,ns] = time.perf_counter() - st
	
		print(n,'\n')
	
	timefig, timeax = plt.subplots()


	timefig.suptitle('Time to test')
	timeax.plot(range(num_samples),mean(elapseds,axis = 0))
	timefig.savefig('stat_figs/timing_'+ str(num_mc_sims) +'_species.png')

	fp, axarrp = plt.subplots(4,2, figsize=(20, 10))
	fpt, axarrpt = plt.subplots(4,2, figsize=(20, 10))
	#fb, ararrb = plt.subplots(4,2, figsize=(20, 10))

	plt.setp(axarrp, xticks=[])
	plt.setp(axarrpt, xticks=[])
	#plt.setp(axarrb, xticks=[])


	fp.suptitle('correlation')
	axarrp[0,0].plot(range(num_samples),mean(num_edges_p,axis = 0))
	axarrp[0,1].plot(range(num_samples),mean(mean_deg_p,axis = 0))
	axarrp[1,1].plot(range(num_samples),mean(p_md_p,axis = 0))
	axarrp[1,0].plot(range(num_samples),mean(p_vd_p,axis = 0))
	axarrp[2,0].plot(range(num_samples),mean(p_max_ev_p,axis = 0))
	axarrp[2,1].plot(range(num_samples),mean(p_min_ev_p,axis = 0))
	axarrp[3,0].plot(range(num_samples),mean(p_num_edges_p,axis = 0))
	axarrp[3,1].axis('off')

	axarrp[0,0].set_title('Number Edges', fontsize=10)
	axarrp[0,1].set_title('Mean Degree', fontsize=10)
	axarrp[1,1].set_title('Mean Degree Probability', fontsize=10)
	axarrp[1,0].set_title('Degree Variance Probability', fontsize=10)
	axarrp[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
	axarrp[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
	axarrp[3,0].set_title('Number Edges Probability', fontsize=10)

	fpt.suptitle('correlation threshholded')
	axarrpt[0,0].plot(range(num_samples),mean(num_edges_pt,axis = 0))
	axarrpt[0,1].plot(range(num_samples),mean(mean_deg_pt,axis = 0))
	axarrpt[1,1].plot(range(num_samples),mean(p_md_pt,axis = 0))
	axarrpt[1,0].plot(range(num_samples),mean(p_vd_pt,axis = 0))
	axarrpt[2,0].plot(range(num_samples),mean(p_max_ev_pt,axis = 0))
	axarrpt[2,1].plot(range(num_samples),mean(p_min_ev_pt,axis = 0))
	axarrpt[3,0].plot(range(num_samples),mean(p_num_edges_pt,axis = 0))
	axarrpt[3,1].axis('off')

	axarrpt[0,0].set_title('Number Edges', fontsize=10)
	axarrpt[0,1].set_title('Mean Degree', fontsize=10)
	axarrpt[1,1].set_title('Mean Degree Probability', fontsize=10)
	axarrpt[1,0].set_title('Degree Variance Probability', fontsize=10)
	axarrpt[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
	axarrpt[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
	axarrpt[3,0].set_title('Number Edges Probability', fontsize=10)

	# fpt.suptitle('binned')
	# axarrb[0,0].plot(range(num_samples),mean(num_edges_b,axis = 0))
	# axarrb[0,1].plot(range(num_samples),mean(mean_deg_b,axis = 0))
	# axarrb[1,1].plot(range(num_samples),mean(p_md_b,axis = 0))
	# axarrb[1,0].plot(range(num_samples),mean(p_vd_b,axis = 0))
	# axarrb[2,0].plot(range(num_samples),mean(p_max_ev_b,axis = 0))
	# axarrb[2,1].plot(range(num_samples),mean(p_min_ev_b,axis = 0))
	# axarrb[3,0].plot(range(num_samples),mean(p_num_edges_b,axis = 0)) 
	# axarrb[3,1].axis('off')
	# 
	# axarrb[0,0].set_title('Number Edges', fontsize=10)
	# axarrb[0,1].set_title('Mean Degree', fontsize=10)
	# axarrb[1,1].set_title('Mean Degree Probability', fontsize=10)
	# axarrb[1,0].set_title('Degree Variance Probability', fontsize=10)
	# axarrb[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
	# axarrb[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
	# axarrb[3,0].set_title('Number Edges Probability', fontsize=10)

	fp.savefig('stat_figs/pears_'+ str(num_mc_sims) +'_species.png')
	fpt.savefig('stat_figs/pears_thresh_'+ str(num_mc_sims) +'_species.png')
	#plt.show()







