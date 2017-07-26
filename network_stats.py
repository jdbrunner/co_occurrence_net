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
		

		
		
		#####run on given networks
		premade = sys.argv[4]
		if premade == 'True':
			premade = True
		else:
			premade = False
		if premade:
			this_array = lvl_abundance_array
			
			#get the adjacency matrix of the network
			adj_name = sys.argv[5]
			type = sys.argv[6]
# 			pears_thr_name = sys.argv[6]
# 			binned_name = sys.argv[7]
			
			adj = pd.read_csv(adj_name, sep = '\t', index_col = 0)#build_network(this_array, 'pearson', thr = False, list_too = False)
# 			pears_thr = pd.read_csv(pears_thr_name, sep = '\t', index_col = 0)#build_network(this_array, 'pearson', thr = True, list_too = False)
# 			binned = pd.read_csv(binned_name, sep = '\t', index_col = 0)#build_network(this_array, 'bins', list_too = False)
# 						
			adj_data = adj.values
# 			pears_thr_data = pears.values
# 			binned_data = binned.values
		
			num_edges = len(adj_data.nonzero()[0])/2
# 			num_edges_pt = len(pears_thr_data.nonzero()[0])/2
# 			num_edges_b = len(binned_data.nonzero()[0])/2
		
			mean_deg = mean(adj_data.sum(axis = 1))
# 			mean_deg_pt = mean(pears_thr_data.sum(axis = 1))
# 			mean_deg_b = mean(binned_data.sum(axis = 1))
			
			this_np_array = this_array.values[:,2:]
			
			fld_where = adj_name.rfind('/')
			fl_nm = adj_name[:fld_where] +'/stats.txt'
			file = open(fl_nm,'a') 
			
			if type == 'pears':
				[p_md, p_vd, p_max_ev, p_min_ev, p_num_edges] = mc_network_stats(this_np_array, adj_data, sims  = 10000)
				file.write('pears: '+i+' '+stype+'\n')
			
			elif type == 'pears_thr':
				[p_md, p_vd, p_max_ev, p_min_ev, p_num_edges] = mc_network_stats(this_np_array, adj_data, thr = True, sims = 10000)
				file.write('pears_thr: '+i+' '+stype+'\n')

			elif type == 'binned':
				[p_md, p_vd, p_max_ev, p_min_ev, p_num_edges] = mc_network_stats(this_np_array, adj_data, bins = True, sims = 10000)
				file.write('binned: '+i+' '+stype+'\n')

			file.write('Number of Edges, Mean Degree, Mean Degree Prob., Degree Variance Prob., Max Eigenvalue Prob., Min Eigenvalue Prob., Number of Edges Prob.\n')
			file.write(str([num_edges, mean_deg, p_md, p_vd, p_max_ev, p_min_ev, p_num_edges])+'\n')
			file.write('All probabilities are probability that random is as high as observed\n\n')

			
			file.close()
			
		build_many = False
		if 	build_many:
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


			for ns in range(1,num_samples):
				print('\r',ns, '/',num_samples - 1, end = '     ')
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
			
				num_edges_p[ns] = len(pears_data.nonzero()[0])/2
				num_edges_pt[ns] = len(pears_thr_data.nonzero()[0])/2
				#num_edges_b[ns] = len(binned_data.nonzero()[0])/2
			
				mean_deg_p[ns] = mean(pears_data.sum(axis = 1))
				mean_deg_pt[ns] = mean(pears_thr_data.sum(axis = 1))
				#mean_deg_b[ns] = mean(binned_data.sum(axis = 1))
			
				num_mc_sims = 1000
				this_np_array = this_array.values[:,2:]
				[p_md_p[ns], p_vd_p[ns], p_max_ev_p[ns], p_min_ev_p[ns],p_num_edges_p[ns]] = mc_network_stats(this_np_array, pears_data, sims = num_mc_sims)
				[p_md_pt[ns], p_vd_pt[ns], p_max_ev_pt[ns], p_min_ev_pt[ns], p_num_edges_pt[ns]] = mc_network_stats(this_np_array, pears_thr_data, thr = True, sims = num_mc_sims)
				#[p_md_b[ns], p_vd_b[ns], p_max_ev_b[ns], p_min_ev_b[ns],p_num_edges_b[ns]] = mc_network_stats(this_np_array, binned_data, bins = True, sims = num_mc_sims)

			fp, axarrp = plt.subplots(4,2, figsize=(20, 10))
			fpt, axarrpt = plt.subplots(4,2, figsize=(20, 10))
			#fb, ararrb = plt.subplots(4,2, figsize=(20, 10))
		
			plt.setp(axarrp, xticks=[])
			plt.setp(axarrpt, xticks=[])
			#plt.setp(axarrb, xticks=[])
		
		
			fp.suptitle('correlation '+stype+' '+i)
			axarrp[0,0].plot(range(num_samples),num_edges_p)
			axarrp[0,1].plot(range(num_samples),mean_deg_p)
			axarrp[1,1].plot(range(num_samples),p_md_p)
			axarrp[1,0].plot(range(num_samples),p_vd_p)
			axarrp[2,0].plot(range(num_samples),p_max_ev_p)
			axarrp[2,1].plot(range(num_samples),p_min_ev_p)
			axarrp[3,0].plot(range(num_samples),p_num_edges_p)
			axarrp[3,1].axis('off')
		
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
			axarrpt[3,1].axis('off')
		
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
	#		axarrb[3,1].axis('off')
		
	# 		axarrb[0,0].set_title('Number Edges', fontsize=10)
	# 		axarrb[0,1].set_title('Mean Degree', fontsize=10)
	# 		axarrb[1,1].set_title('Mean Degree Probability', fontsize=10)
	# 		axarrb[1,0].set_title('Degree Variance Probability', fontsize=10)
	# 		axarrb[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
	# 		axarrb[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
	# 		axarrb[3,0].set_title('Number Edges Probability', fontsize=10)

			fp.savefig('stat_figs/pears_'+ str(num_mc_sims) +'_' +stype +'_'+ i+ '.png')
			fpt.savefig('stat_figs/pears_thresh_'+ str(num_mc_sims) +'_' +stype +'_'+  i+  '.png')
			#fb.savefig('stat_figs/binned_'+ str(num_mc_sims) + '_' +stype +'_'+  i+ '.png')
		
			#plt.show()







