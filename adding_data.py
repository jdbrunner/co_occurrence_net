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
from joblib import Parallel, delayed, cpu_count


if __name__ == '__main__':

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

	monte = False#monte carlo is slow!
	pltrows = 1
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

	#Combine samples of the same kind (that is, the same part of the body) so that we can 
	#color according to where abundance of a genome is highest.
	diff_samps_types = unique([name[:-10] for name in array(abundance_array_full.keys())[2:]])
	for sty in range(len(diff_samps_types)):
		if 'female' in diff_samps_types[sty]:
			diff_samps_types[sty] = diff_samps_types[sty][:-7]
		elif 'male' in diff_samps_types[sty]:
			diff_samps_types[sty] = diff_samps_types[sty][:-5]
	diff_samps_types = unique(diff_samps_types)

	ab_arrays = [abundance_array_full]	
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

	fld_save = 'merged2'
		
	stype_ns = ['Full_'] + [spty + '_' for spty in diff_samps_types]
	for i in level:
		for ii in range(len(ab_arrays)):
			if sample_types:
				stype = stype_ns[ii]
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
			
			num_to_make = len(arange(10,num_samples,10))
			
			num_edges_p = zeros(num_to_make)
			mean_deg_p = zeros(num_to_make)
			if monte:
				p_md_p = zeros(num_to_make)
				p_vd_p = zeros(num_to_make)
				p_max_ev_p = zeros(num_to_make)
				p_min_ev_p = zeros(num_to_make)
				p_num_edges_p = zeros(num_to_make)
	
			num_edges_pt = zeros(num_to_make)
			mean_deg_pt = zeros(num_to_make)
			if monte:
				p_md_pt = zeros(num_to_make)
				p_vd_pt = zeros(num_to_make)
				p_max_ev_pt = zeros(num_to_make)
				p_min_ev_pt = zeros(num_to_make)
				p_num_edges_pt = zeros(num_to_make)
	
	# 		num_edges_b = zeros(num_samples)
	# 		mean_deg_b = zeros(num_samples)
	#		if monte:
		# 		p_md_b = zeros(num_samples)
		# 		p_vd_b = zeros(num_samples)
		# 		p_max_ev_b = zeros(num_samples)
		# 		p_min_ev_b = zeros(num_samples)	
		#		p_num_edges_b = zeros(num_samples)

			ar_count = 0
			for ns in arange(10,num_samples,10):
				print('\r',ns, '/',num_samples, end = '     ')
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
		
				num_edges_p[ar_count] = len(pears_data.nonzero()[0])/2
				num_edges_pt[ar_count] = len(pears_thr_data.nonzero()[0])/2
				#num_edges_b[ar_count] = len(binned_data.nonzero()[0])/2
		
				mean_deg_p[ar_count] = mean(pears_data.sum(axis = 1))
				mean_deg_pt[ar_count] = mean(pears_thr_data.sum(axis = 1))
				#mean_deg_b[ar_count] = mean(binned_data.sum(axis = 1))
				
				
				if monte:		
					num_mc_sims = 1000
					this_np_array = this_array.values[:,2:]
					[p_md_p[ar_count], p_vd_p[nar_counts], p_max_ev_p[ar_count], p_min_ev_p[ar_count],p_num_edges_p[ar_count]] = mc_network_stats(this_np_array, pears_data, sims = num_mc_sims)
					[p_md_pt[ar_count], p_vd_pt[ar_count], p_max_ev_pt[ar_count], p_min_ev_pt[ar_count], p_num_edges_pt[ar_count]] = mc_network_stats(this_np_array, pears_thr_data, thr = True, sims = num_mc_sims)
					#[p_md_b[ns], p_vd_b[ns], p_max_ev_b[ns], p_min_ev_b[ns],p_num_edges_b[ns]] = mc_network_stats(this_np_array, binned_data, bins = True, sims = num_mc_sims)
				
				ar_count = ar_count + 1
				
			if monte:
				pltrows = 4
			fp, axarrp = plt.subplots(pltrows,2, figsize=(20, 10))
			fpt, axarrpt = plt.subplots(pltrows,2, figsize=(20, 10))
			#fb, ararrb = plt.subplots(pltrows,2, figsize=(20, 10))
	
			plt.setp(axarrp, xticks=[])
			plt.setp(axarrpt, xticks=[])
			#plt.setp(axarrb, xticks=[])
	
	
			fp.suptitle('correlation '+stype+' '+i)
			axarrp[0].plot(arange(10,num_samples,10),num_edges_p)
			axarrp[1].plot(arange(10,num_samples,10),mean_deg_p)
			axarrp[0].set_title('Number Edges', fontsize=10)
			axarrp[1].set_title('Mean Degree', fontsize=10)
			if monte:
				axarrp[0,0].plot(arange(10,num_samples,10),num_edges_p)
				axarrp[0,1].plot(arange(10,num_samples,10),mean_deg_p)
				axarrp[0,0].set_title('Number Edges', fontsize=10)
				axarrp[0,1].set_title('Mean Degree', fontsize=10)
				axarrp[1,1].plot(arange(10,num_samples,10),p_md_p)
				axarrp[1,0].plot(arange(10,num_samples,10),p_vd_p)
				axarrp[2,0].plot(arange(10,num_samples,10),p_max_ev_p)
				axarrp[2,1].plot(arange(10,num_samples,10),p_min_ev_p)
				axarrp[3,0].plot(arange(10,num_samples,10),p_num_edges_p)
				axarrp[3,1].axis('off')
	

				axarrp[1,1].set_title('Mean Degree Probability', fontsize=10)
				axarrp[1,0].set_title('Degree Variance Probability', fontsize=10)
				axarrp[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
				axarrp[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
				axarrp[3,0].set_title('Number Edges Probability', fontsize=10)
	
			fpt.suptitle('correlation threshholded '+stype+' '+i)
			axarrpt[0].plot(arange(10,num_samples,10),num_edges_pt)
			axarrpt[1].plot(arange(10,num_samples,10),mean_deg_pt)
			axarrpt[0].set_title('Number Edges', fontsize=10)
			axarrpt[1].set_title('Mean Degree', fontsize=10)
			if monte:
				axarrpt[0,0].plot(arange(10,num_samples,10),num_edges_pt)
				axarrpt[0,1].plot(arange(10,num_samples,10),mean_deg_pt)
				axarrpt[0,0].set_title('Number Edges', fontsize=10)
				axarrpt[0,1].set_title('Mean Degree', fontsize=10)
				axarrpt[1,1].plot(arange(10,num_samples,10),p_md_pt)
				axarrpt[1,0].plot(arange(10,num_samples,10),p_vd_pt)
				axarrpt[2,0].plot(arange(10,num_samples,10),p_max_ev_pt)
				axarrpt[2,1].plot(arange(10,num_samples,10),p_min_ev_pt)
				axarrpt[3,0].plot(arange(10,num_samples,10),p_num_edges_pt)
				axarrpt[3,1].axis('off')
	

				axarrpt[1,1].set_title('Mean Degree Probability', fontsize=10)
				axarrpt[1,0].set_title('Degree Variance Probability', fontsize=10)
				axarrpt[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
				axarrpt[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
				axarrpt[3,0].set_title('Number Edges Probability', fontsize=10)
	
	# 		fpt.suptitle('binned '+stype+' '+i)
	# 		axarrb[0].plot(range(num_samples),num_edges_b)
	# 		axarrb[1].plot(range(num_samples),mean_deg_b)
		# 	axarrb[0].set_title('Number Edges', fontsize=10)
	# 		axarrb[1].set_title('Mean Degree', fontsize=10)
# 			if monte:
		# 		axarrb[0,0].plot(range(num_samples),num_edges_b)
		# 		axarrb[0,1].plot(range(num_samples),mean_deg_b)
			# 	axarrb[0,0].set_title('Number Edges', fontsize=10)
		# 		axarrb[0,1].set_title('Mean Degree', fontsize=10)
		# 		axarrb[1,1].plot(range(num_samples),p_md_b)
		# 		axarrb[1,0].plot(range(num_samples),p_vd_b)
		# 		axarrb[2,0].plot(range(num_samples),p_max_ev_b)
		# 		axarrb[2,1].plot(range(num_samples),p_min_ev_b)
		# 		axarrb[3,0].plot(range(num_samples),p_num_edges_b) 
		#		axarrb[3,1].axis('off')
	

		# 		axarrb[1,1].set_title('Mean Degree Probability', fontsize=10)
		# 		axarrb[1,0].set_title('Degree Variance Probability', fontsize=10)
		# 		axarrb[2,0].set_title('Max Eigenvalue Probability', fontsize=10)
		# 		axarrb[2,1].set_title('Min Eigenvalue Probability', fontsize=10)
		# 		axarrb[3,0].set_title('Number Edges Probability', fontsize=10)
			
			if monte:
				fp.savefig('stat_figs/'+ fld_save + '/pears_'+ str(num_mc_sims) +'_' +stype +'_'+ i+ '.png')
				fpt.savefig('stat_figs/'+ fld_save + '/pears_thresh_'+ str(num_mc_sims) +'_' +stype +'_'+  i+  '.png')
				#fb.savefig('stat_figs/'+ fld_save + '/binned_'+ str(num_mc_sims) + '_' +stype +'_'+  i+ '.png')
			
			else:
				fp.savefig('stat_figs/'+ fld_save + '/pears_' +stype +'_'+ i+ '.png')
				fpt.savefig('stat_figs/'+ fld_save + '/pears_thresh_'+stype +'_'+  i+  '.png')
				#fb.savefig('stat_figs/'+ fld_save + '/binned_'+stype +'_'+  i+ '.png')
			
			plt.close(fp)
			plt.close(fpt)
			
			#plt.show()







