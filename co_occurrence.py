#########################################################################################
#																						#		
#				Co-Occurrence Network Builder											#
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
import itertools as iter
from scipy.special import binom as choose
from scipy.stats import binom
from numpy.random import binomial as bino
import re
import time
import numpy.ma as mask
from co_occ_funs import *

if __name__ == '__main__':


	#####Variables
	##User input, should be entered as an argument
	csv_name = sys.argv[1]#the data
	level = sys.argv[2]#which taxonomic lvls do you care about
	net_name = sys.argv[3]#folder you want it saved in
	sample_types = sys.argv[4]#should it separate by sample type
	if sample_types == 'True':
		sample_types = True
	else:
		sample_types = False
	
	numhld = int(sys.argv[5])#hold out some columns to test on
	if numhld != 0:
		hldouts = True
	else:
		hldouts = False

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
	#print(abundance_array)
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



	#Combine samples of the same kind (that is, the same part of the body) so that we can 
	#color according to where abundance of a genome is highest.
	diff_samps_types = unique([name[:-10] for name in array(abundance_array_full.keys())[2:]])
	for sty in range(len(diff_samps_types)):
		if 'female' in diff_samps_types[sty]:
			diff_samps_types[sty] = diff_samps_types[sty][:-7]
		elif 'male' in diff_samps_types[sty]:
			diff_samps_types[sty] = diff_samps_types[sty][:-5]
	diff_samps_types = unique(diff_samps_types)
	data_samps_types = transpose([sum([abundance_array_full[samp] for samp in abundance_array_full.keys()[2:]
								 if smp_type in samp],axis = 0) for smp_type in diff_samps_types])

	samp_type_abund = pd.DataFrame(data_samps_types, columns = diff_samps_types)

	samp_type_norm = True
	if samp_type_norm:
		for col in samp_type_abund.columns:
			tot = sum(samp_type_abund[col])
			samp_type_abund[col] = samp_type_abund[col]/tot


	samp_type_abund['TAXA'] = abundance_array_full['TAXA']

	the_levels = abundance_array_full['LEVEL'].values
	the_taxa = abundance_array_full['TAXA'].values


	gender = True
	other_meta = False
	
	ab_arrays = [abundance_array_full]	
	#gender and sample location are in the sample titles. For other meta, we need a metadata file
	if sample_types:
		for samp in diff_samps_types:
			ab_samp = pd.DataFrame(the_levels, columns = ['LEVEL'], index = abundance_array_full.index)
			ab_samp['TAXA'] = the_taxa
			selecter = abundance_array_full.columns.map(lambda x: bool(re.search(samp,x)))
			ab_samp[abundance_array_full.columns[where(selecter)]] = abundance_array_full[abundance_array_full.columns[where(selecter)]]
			ab_arrays = ab_arrays + [ab_samp]
	elif gender:
		for gen in ['male','female']:
			ab_gen = pd.DataFrame(the_levels, columns = ['LEVEL'], index = abundance_array_full.index)
			ab_gen['TAXA'] = the_taxa
			selecter = abundance_array_full.columns.map(lambda x: bool(re.search(gen,x)))
			ab_gen[abundance_array_full.columns[where(selecter)]] = abundance_array_full[abundance_array_full.columns[where(selecter)]]
			ab_arrays = ab_arrays + [ab_gen]
	elif other_meta:
		meta_file = sys.argv[6]
		meta_column = sys.argv[7]
		meta_arr = pd.read_csv(meta_file, sep = ' ', index_col= 0)
		meta_col_dat = meta_arr.meta_column
		meta_groups = unique(meta_col_dat.values)
		if len(meta_groups) > 0.9*len(meta_arr):
			#if the category is numerical (BMI for example), we need to bin it
			numbins = 10
			meta_groups = linspace(min(meta_col_dat), max(meta_col_dat), numbins)
			x = subtract.outer(meta_col_dat, meta_groups) #row j is meta_col_dat[j] - meta_groups
			y = argmin(abs(x),axis = 1) #entry j is entry number of meta_groups closest to meta_col_dat[j]
			catdat = pd.DataFrame(meta_groups[y], columns = [c], index = meta_arr.index)
		else:
			catdat = meta_col_dat
		for cl in meta_groups:
			id_list = meta_arr.index[meta_col_dat == cl].values
			ab_met = pd.DataFrame(the_levels, columns = ['LEVEL'], index = abundance_array_full.index)
			ab_met['TAXA'] = the_taxa
			selecter = where([any([id in colnm for id in id_list]) for colnm in abundance_array_full.columns])
			ab_met[abundance_array_full.columns[selecter]] = abundance_array_full[abundance_array_full.columns[selecter]]
			ab_arrays = ab_arrays + [ab_met]


	if hldouts:
		L = len(abundance_array_full.columns) - 2
		hld = randint(L, size = numhld) + 2
		abundance_array_full.drop(abundance_array_full.columns[hld],axis = 1, inplace = True)

	save = True
	stype_ns = ['Full_'] + [spty + '_' for spty in diff_samps_types]
	gtype = ['Full_','Male_','Female_']
	for i in level:
		print(i)
		for ii in range(len(ab_arrays)):
			if sample_types:
				stype = stype_ns[ii]
			elif gender:
				
			else:
				stype = ''
			abundance_array = ab_arrays[ii]
			print(stype)

		
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
				
			[adjacency_frames_bins, source_target_frames_bins] = build_network(lvl_abundance_array, 'bins')
			[adjacency_frames_pear, source_target_frames_pear] = build_network(lvl_abundance_array, 'pearson')
			[adjacency_frames_pear_thr, source_target_frames_pear_thr] = build_network(lvl_abundance_array, 'pearson', thr = True)

		
			[source_target_frames_bins, node_data_bins] = make_meta(source_target_frames_bins, samp_type_abund, lvl_abundance_array)
			[source_target_frames_pear, node_data_pear] = make_meta(source_target_frames_pear, samp_type_abund, lvl_abundance_array)
			[source_target_frames_pear_thr, node_data_pear_thr] = make_meta(source_target_frames_pear_thr, samp_type_abund, lvl_abundance_array)


			if save:
					
					
				flname1 = net_name+'/bins/'+i + '_' + stype+'_adj.tsv'
				flname2 = net_name+'/bins/'+i + '_' + stype+'_list.tsv'
				flname5 = net_name+'/bins/'+i + '_' + stype+'_node_data.tsv'
			
				flname3t = net_name+'/pears/'+i + '_' + stype+'thr_adj.tsv'
				flname4t = net_name+'/pears/'+i + '_' + stype+'thr_list.tsv'
				flname6t = net_name+'/pears/'+i + '_' + stype+'thr_node_data.tsv'
		
				flname3 = net_name+'/pears/'+i + '_' + stype+'cor_adj.tsv'
				flname4 = net_name+'/pears/'+i + '_' + stype+'cor_list.tsv'
				flname6 = net_name+'/pears/'+i + '_' + stype+'cor_node_data.tsv'
				
				if hldouts:
					holddf = pd.DataFrame(hld, columns = ['Holdouts'])
					hldfile1 = flname1[:-8] + '_held.tsv' 
					hldfile2 = flname3t[:-8] + '_held.tsv'
					hldfile3 = flname3[:-8] + '_held.tsv'
					holddf.to_csv(hldfile1, sep = '\t')
					holddf.to_csv(hldfile2, sep = '\t')
					holddf.to_csv(hldfile3, sep = '\t')
					
				
			
				adjacency_frames_bins.to_csv(flname1, sep = '\t')
				source_target_frames_bins.to_csv(flname2, sep = '\t', index = False)
				node_data_bins.to_csv(flname5, sep = '\t')

				adjacency_frames_pear.to_csv(flname3, sep = '\t')
				source_target_frames_pear.to_csv(flname4, sep = '\t', index = False)
				node_data_pear.to_csv(flname6, sep = '\t')
			
				adjacency_frames_pear_thr.to_csv(flname3t, sep = '\t')
				source_target_frames_pear_thr.to_csv(flname4t, sep = '\t', index = False)
				node_data_pear_thr.to_csv(flname6t, sep = '\t')

















