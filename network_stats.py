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


#import networkx as nx

#####Variables
##User input, should be entered as an argument
csv_name = sys.argv[1]#the data
level = sys.argv[2]
sample_types = sys.argv[3]#should it separate by sample type
if sample_types == 'True':
	sample_types = True
else:
	sample_types = False


##I have to parse the level input I guess...
if level != 'all':
	lvl_list = []
	if level[0] == '[' and level[-1] == ']':
		commas = argwhere([let == ',' for let in level])
		word_start = 1
		for com in commas:
			lvl_list = lvl_list + [level[word_start:com[0]]]
			word_start = com[0]+1
		lvl_list = lvl_list+[level[word_start:-1]]
	else:
		lvl_list  = [level]
	level = lvl_list


######Import the abundance matrix
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 
abundance_array_full = pd.read_csv(csv_name, sep = ' ')
#Let's remove taxa that aren't present at all.
for ind in abundance_array_full.index:
	if max(abundance_array_full.loc[ind][2:]) == 0:
		abundance_array_full = abundance_array_full.drop(ind)

abundance_array_full.index = range(len(abundance_array_full))		

#We want to make a network at a specific level, or set of levels
if level == 'all':
	level = array(abundance_array_full['LEVEL'].values)
	level = array(unique(level))
else:
	level = array(level)

#need to check that the levels given as argument actually exist
levels_allowed = array(abundance_array_full['LEVEL'].values)
levels_allowed = array(unique(levels_allowed))

for i in level:
	if not i in levels_allowed:
		print(levels_allowed)
		sys.exit('Bad level argument. Choose a level from the list above or type all')


the_levels = abundance_array_full['LEVEL'].values
the_taxa = abundance_array_full['TAXA'].values

ab_arrays = []	
if sample_types:
	for samp in diff_samps_types:
		slist = []
		ab_samp = pd.DataFrame(the_levels, columns = ['LEVEL'], index = abundance_array_full.index)
		ab_samp['TAXA'] = the_taxa
		selecter = abundance_array_full.columns.map(lambda x: bool(re.search(samp,x)))
		ab_samp[abundance_array_full.columns[where(selecter)]] = abundance_array_full[abundance_array_full.columns[where(selecter)]]
		ab_arrays = ab_arrays + [ab_samp]
else:
	ab_arrays = [abundance_array_full]


for i in level:
	for ii in range(len(ab_arrays)):
		if sample_types:
			stype = diff_samps_types[ii]
		else:
			stype = ''
		abundance_array = ab_arrays[ii]
		#start by getting the indices of the members of the level
		in_level = where(abundance_array['LEVEL'] == i)[0]
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
		

		
		
		#####TESTING
		testing = False
		if testing:
			this_array = lvl_abundance_array
			pears = build_network(this_array, 'pearson', thr = False, list_too = False)
			pears_thr = build_network(this_array, 'pearson', thr = True, list_too = False)
			binned = build_network(this_array, 'bins', list_too = False)
			
			print('built')
			
			pears_data = pears.values
			pears_thr_data = pears.values
			binned_data = binned.values
		
			num_edges_p = len(pears_data.nonzero()[0])/2
			num_edges_pt = len(pears_thr_data.nonzero()[0])/2
			num_edges_b = len(binned_data.nonzero()[0])/2
		
			mean_deg_p = mean(pears_data.sum(axis = 1))
			mean_deg_pt = mean(pears_thr_data.sum(axis = 1))
			mean_deg_b = mean(binned_data.sum(axis = 1))
			
			this_np_array = this_array.values[:,2:]
			[p_md_p, p_vd_p, p_max_ev_p, p_min_ev_p] = mc_network_stats(this_np_array, pears_data)
			print('pears')
			print([num_edges_p, mean_deg_p, p_md_p, p_vd_p, p_max_ev_p, p_min_ev_p])
			
			[p_md_pt, p_vd_pt, p_max_ev_pt, p_min_ev_pt] = mc_network_stats(this_np_array, pears_thr_data, thr = True)
			print('pears_thr')
			print([num_edges_pt, mean_deg_pt, p_md_pt, p_vd_pt, p_max_ev_pt, p_min_ev_pt])
			
			[p_md_b, p_vd_b, p_max_ev_b, p_min_ev_b] = mc_network_stats(this_np_array, binned_data, bins = True)
			print('binned')
			print([num_edges_b, mean_deg_b, p_md_b, p_vd_b, p_max_ev_b, p_min_ev_b])
	
			sys.exit()
			
		num_edges_p = zeros(num_samples)
		mean_deg_p = zeros(num_samples)
		p_md_p = zeros(num_samples)
		p_vd_p = zeros(num_samples)
		p_max_ev_p = zeros(num_samples)
		p_min_ev_p = zeros(num_samples)
		p_num_edges_p = zeros(num_samples)
		
		num_edges_pt = zeros(num_samples)
		mean_deg_pt = zeros(num_samples)
		p_md_pt = zeros(num_samples)
		p_vd_pt = zeros(num_samples)
		p_max_ev_pt = zeros(num_samples)
		p_min_ev_pt = zeros(num_samples)
		p_num_edges_pt = zeros(num_samples)
		
# 		num_edges_b = zeros(num_samples)
# 		mean_deg_b = zeros(num_samples)
# 		p_md_b = zeros(num_samples)
# 		p_vd_b = zeros(num_samples)
# 		p_max_ev_b = zeros(num_samples)
# 		p_min_ev_b = zeros(num_samples)	
#		p_num_edges_b = zeros(num_samples)


		for ns in range(1,10):#num_samples):
			#print('\r',ns, '/',num_samples - 1, end = '     ')
			chs = array(sample(range(num_samples), ns)) + 2
			this_array = lvl_abundance_array.loc[:,lvl_abundance_array.columns[append([0,1],chs)]]
			print(this_array.columns)
			#build a network with a subset of the data
			pears = build_network(this_array, 'pearson', thr = False, list_too = False)
			pears_thr = build_network(this_array, 'pearson', thr = True, list_too = False)
			#binned = build_network(this_array, 'bins', list_too = False)
			
			pears_data = pears.values
			pears_thr_data = pears_thr.values
			#binned_data = binned.values
			
			num_edges_p[ns] = len(pears_data.nonzero()[0])/2
			num_edges_pt[ns] = len(pears_thr_data.nonzero()[0])/2
			#num_edges_b[ns] = len(binned_data.nonzero()[0])/2
			
			mean_deg_p[ns] = mean(pears_data.sum(axis = 1))
			mean_deg_pt[ns] = mean(pears_thr_data.sum(axis = 1))
			#mean_deg_b[ns] = mean(binned_data.sum(axis = 1))
			
			
			this_np_array = this_array.values[:,2:]
			[p_md_p[ns], p_vd_p[ns], p_max_ev_p[ns], p_min_ev_p[ns],p_num_edges_p[ns]] = mc_network_stats(this_np_array, pears_data, sims = 1000)
			[p_md_pt[ns], p_vd_pt[ns], p_max_ev_pt[ns], p_min_ev_pt[ns], p_num_edges_pt[ns]] = mc_network_stats(this_np_array, pears_thr_data, thr = True, sims = 1000)
			#[p_md_b[ns], p_vd_b[ns], p_max_ev_b[ns], p_min_ev_b[ns],p_num_edges_b[ns]] = mc_network_stats(this_np_array, binned_data, bins = True, sims = 1000)

		fp, axarrp = plt.subplots(4,2)
		fpt, axarrpt = plt.subplots(4,2)
		#fb, ararrb = plt.subplots(4,2)
		
		fp.suptitle('correlation '+stype+' '+i)
		axarrp[0,0].plot(range(num_samples),num_edges_p)
		axarrp[0,1].plot(range(num_samples),mean_deg_p)
		axarrp[1,1].plot(range(num_samples),p_md_p)
		axarrp[1,0].plot(range(num_samples),p_vd_p)
		axarrp[2,0].plot(range(num_samples),p_max_ev_p)
		axarrp[2,1].plot(range(num_samples),p_min_ev_p)
		axarrp[3,0].plot(range(num_samples),p_num_edges_p)
		
		axarrp[0,0].set_title('Number Edges', fontsize=10)
		axarrp[0,1].set_title('Mean Degree', fontsize=10)
		axarrp[1,1].set_title('Mean Degree Probability', fontsize=10)
		axarrp[1,0].set_title('Degree Variance Probability', fontsize=10)
		axarrp[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
		axarrp[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
		axarrp[3,0].set_title('Number Edges Probability', fontsize=10)
		
		fpt.suptitle('correlation threshholded '+stype+' '+i)
		axarrpt[0,0].plot(range(num_samples),num_edges_pt)
		axarrpt[0,1].plot(range(num_samples),mean_deg_pt)
		axarrpt[1,1].plot(range(num_samples),p_md_pt)
		axarrpt[1,0].plot(range(num_samples),p_vd_pt)
		axarrpt[2,0].plot(range(num_samples),p_max_ev_pt)
		axarrpt[2,1].plot(range(num_samples),p_min_ev_pt)
		axarrpt[3,0].plot(range(num_samples),p_num_edges_pt)
		
		axarrpt[0,0].set_title('Number Edges', fontsize=10)
		axarrpt[0,1].set_title('Mean Degree', fontsize=10)
		axarrpt[1,1].set_title('Mean Degree Probability', fontsize=10)
		axarrpt[1,0].set_title('Degree Variance Probability', fontsize=10)
		axarrpt[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
		axarrpt[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
		axarrpt[3,0].set_title('Number Edges Probability', fontsize=10)
		
# 		fpt.suptitle('binned '+stype+' '+i)
# 		axarrb[0,0].plot(range(num_samples),num_edges_b)
# 		axarrb[0,1].plot(range(num_samples),mean_deg_b)
# 		axarrb[1,1].plot(range(num_samples),p_md_b)
# 		axarrb[1,0].plot(range(num_samples),p_vd_b)
# 		axarrb[2,0].plot(range(num_samples),p_max_ev_b)
# 		axarrb[2,1].plot(range(num_samples),p_min_ev_b)
# 		axarrb[3,0].plot(range(num_samples),p_num_edges_b)
# 		
# 		axarrb[0,0].set_title('Number Edges', fontsize=10)
# 		axarrb[0,1].set_title('Mean Degree', fontsize=10)
# 		axarrb[1,1].set_title('Mean Degree Probability', fontsize=10)
# 		axarrb[1,0].set_title('Degree Variance Probability', fontsize=10)
# 		axarrb[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
# 		axarrb[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
# 		axarrb[3,0].set_title('Number Edges Probability', fontsize=10)

		fp.savefig('stat_figs/pears.png')
		fpt.savefig('stat_figs/pears_thresh.png')
		#fb.savefig('stat_figs/binned.png')
		
		#plt.show()







