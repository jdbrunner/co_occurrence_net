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

hldouts = sys.argv[5]#hold out some columns to test on
if hldouts == 'True':
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


if hldouts:
	L = len(abundance_array_full.columns) - 2
	hld = randint(L, size = 4) + 2
	abundance_array_full.drop(abundance_array_full.columns[hld],axis = 1, inplace = True)

save = True
for i in level:
	print(i)
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
				
		[adjacency_frames_bins, source_target_frames_bins] = build_network(lvl_abundance_array, 'bins')
		[adjacency_frames_pear, source_target_frames_pear] = build_network(lvl_abundance_array, 'pearson')
		[adjacency_frames_pear_thr, source_target_frames_pear_thr] = build_network(lvl_abundance_array, 'pearson', thr = True)

 		
		[source_target_frames_bins, node_data_bins] = make_meta(source_target_frames_bins, samp_type_abund, lvl_abundance_array)
		[source_target_frames_pear, node_data_pear] = make_meta(source_target_frames_pear, samp_type_abund, lvl_abundance_array)
		[source_target_frames_pear_thr, node_data_pear_thr] = make_meta(source_target_frames_pear_thr, samp_type_abund, lvl_abundance_array)


		if save:
			if hldouts:
				held = '_'.join(map(str,hld))
				flname1 = net_name+'/bins/'+i + '_' + stype+'_'+held+'_adj.tsv'
				flname2 = net_name+'/bins/'+i + '_' + stype+'_'+held+'_list.tsv'
				flname5 = net_name+'/bins/'+i + '_' + stype+'_'+held+'_node_data.tsv'
				
				flname3t = net_name+'/pears/'+i + '_' + stype+'_'+held+'thr_adj.tsv'
				flname4t = net_name+'/pears/'+i + '_' + stype+'_'+held+'thr_list.tsv'
				flname6t = net_name+'/pears/'+i + '_' + stype+'_'+held+'thr_node_data.tsv'
			
				flname3 = net_name+'/pears/'+i + '_' + stype+'_'+held+'cor_adj.tsv'
				flname4 = net_name+'/pears/'+i + '_' + stype+'_'+held+'cor_list.tsv'
				flname6 = net_name+'/pears/'+i + '_' + stype+'_'+held+'cor_node_data.tsv'
			else:
				flname1 = net_name+'/bins/'+i + '_' + stype+'_adj.tsv'
				flname2 = net_name+'/bins/'+i + '_' + stype+'_list.tsv'
				flname5 = net_name+'/bins/'+i + '_' + stype+'_node_data.tsv'
				
				flname3t = net_name+'/pears/'+i + '_' + stype+'thr_adj.tsv'
				flname4t = net_name+'/pears/'+i + '_' + stype+'thr_list.tsv'
				flname6t = net_name+'/pears/'+i + '_' + stype+'thr_node_data.tsv'
			
				flname3 = net_name+'/pears/'+i + '_' + stype+'cor_adj.tsv'
				flname4 = net_name+'/pears/'+i + '_' + stype+'cor_list.tsv'
				flname6 = net_name+'/pears/'+i + '_' + stype+'cor_node_data.tsv'
				
			
			adjacency_frames_bins.to_csv(flname1, sep = '\t')
			source_target_frames_bins.to_csv(flname2, sep = '\t')
			node_data_bins.to_csv(flname5, sep = '\t')

			adjacency_frames_pear.to_csv(flname3, sep = '\t')
			source_target_frames_pear.to_csv(flname4, sep = '\t')
			node_data_pear.to_csv(flname6, sep = '\t')
			
			adjacency_frames_pear_thr.to_csv(flname3t, sep = '\t')
			source_target_frames_pear_thr.to_csv(flname4t, sep = '\t')
			node_data_pear_thr.to_csv(flname6t, sep = '\t')

















