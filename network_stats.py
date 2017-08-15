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
		

		
		
			#####run on given networks
			this_array = lvl_abundance_array
		
			#get the adjacency matrix of the network
			adj_name = sys.argv[4]
			type = sys.argv[5]
			
			go_on = True
			#only compare type to type
			if sample_types:
				if not(stype in adj_name):
					go_on = False


			if go_on:
				adj = pd.read_csv(adj_name, sep = '\t', index_col = 0)#build_network(this_array, 'pearson', thr = False, list_too = False)

				adj_data = adj.values

	
				num_edges = len(adj_data.nonzero()[0])/2

				if num_edges > 0:
				
	
					mean_deg = mean(adj_data.sum(axis = 1))

					this_np_array = this_array.values[:,2:]
		
					fld_where = adj_name.rfind('/')
					fl_nm = adj_name[:fld_where] +'/stats.txt'
					file = open(fl_nm,'a') 
			
	# 				print(adj_name)
	# 				print(type)
		
					if type == 'pears':
						[p_md, p_vd, p_max_ev, p_min_ev, p_num_edges] = mc_network_stats(this_np_array, adj_data, sims  = 100)
						file.write('pears: '+i+' '+stype+'\n')
		
					elif type == 'pears_thr':
						[p_md, p_vd, p_max_ev, p_min_ev, p_num_edges] = mc_network_stats(this_np_array, adj_data, thr = True, sims = 100)
						file.write('pears_thr: '+i+' '+stype+'\n')

					elif type == 'binned':
						[p_md, p_vd, p_max_ev, p_min_ev, p_num_edges] = mc_network_stats(this_np_array, adj_data, bins = True, sims = 100)
						file.write('binned: '+i+' '+stype+'\n')

					file.write('Number of Edges, Mean Degree, Mean Degree Prob., Degree Variance Prob., Max Eigenvalue Prob., Min Eigenvalue Prob., Number of Edges Prob.\n')
					file.write(str([num_edges, mean_deg, p_md, p_vd, p_max_ev, p_min_ev, p_num_edges])+'\n')
					file.write('All probabilities are probability that random is as high as observed\n\n')

		
					file.close()
		
		






