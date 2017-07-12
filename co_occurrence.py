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
csv_name = sys.argv[1]
level = sys.argv[2]
net_name = sys.argv[3]
sample_types = sys.argv[4]
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


save = True
old_way = True
new_way = True
# adjacency_frames_bins = dict()
# source_target_frames_bins = dict()
# adjacency_frames_pear = dict()
# source_target_frames_pear = dict()
#colors = dict()
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
		lvl_abundance_array = lvl_abundance_array.drop('LEVEL',axis = 1)
		#create np array of abundance values (dropping labels) within the level
		ab_np_array1 = lvl_abundance_array.values[:,1:]
		not_seen = where(sum(ab_np_array1,1) == 0)[0]
		n_seen_ind = lvl_abundance_array.index[not_seen]
		in_level = delete(in_level,not_seen)
		#prune off unseen organisms
		lvl_abundance_array = lvl_abundance_array.drop(n_seen_ind)
		ab_np_array = lvl_abundance_array.values[:,1:].astype(float)

		
		if old_way:
	#Create numpy array of co-occurrence fractions. This is the incidence matrix for our weighted
	#graph.
		#t1 = time.time()
			occur_probs = occ_probs(lvl_abundance_array,0.05,5)
			adj_matrix_pre_bins = asarray([[0 if x == y else both_same_bin(ab_np_array[x],ab_np_array[y], 0.05, 5) 
										for x in range(len(ab_np_array))] for y in range(len(ab_np_array))])
			
			degs_b = sum(adj_matrix_pre_bins,1)
			unconnected_b = where(degs_b == 0)
			adj_matrix_pre_bins = delete(adj_matrix_pre_bins,unconnected_b,0)
			adj_matrix_pre_bins = delete(adj_matrix_pre_bins,unconnected_b,1)
			ab_np_array_bins = delete(ab_np_array,unconnected_b,0)
			if ab_np_array_bins.shape[0] != 0:
				in_level_b = delete(in_level,unconnected_b)
			
				stringency_b = 0.05
				N = len(ab_np_array[0])
				tot_seen_b = [sum([abd != 0 for abd in row]) for row in ab_np_array]
			
				adj_size_b = adj_matrix_pre_bins.shape
				adj_matrix_bins = zeros(adj_size_b)
				for k in range(adj_size_b[0]):
					for j in range(adj_size_b[1]):
						if adj_matrix_pre_bins[k,j] != 0:
							p = approx_rand_prob(occur_probs,adj_matrix_pre_bins[k,j],k,j)
							if p <= stringency_b:
								adj_matrix_bins[k,j] = adj_matrix_pre_bins[k,j] 
			
				degs2_b = sum(adj_matrix_bins,1)
				unconnected2_b = where(degs2_b == 0)
				adj_matrix_bins = delete(adj_matrix_bins,unconnected2_b,0)
				adj_matrix_bins = delete(adj_matrix_bins,unconnected2_b,1)
				tot_seen_2_b = delete(tot_seen_b,unconnected2_b)
				
				in_level_2_b = delete(in_level_b,unconnected2_b)				
							
				adjacency_frames_bins = pd.DataFrame(adj_matrix_bins, index = abundance_array['TAXA'][in_level_2_b], columns = abundance_array['TAXA'][in_level_2_b])


				#create a list of edges - source in one column, target in the other. This is an alternative to the adjacency matrix
				#that might make it easier to load in to cytoscape.
				num_edge_b = count_nonzero(adj_matrix_bins)
				source_target_data_bins = []#empty([num_edge,3], dtype = 'string')	
				for l in range(len(adj_matrix_bins)):
					for k in range(l+1):
							if adj_matrix_bins[l,k] != 0:
								s_type1 = color_picker(samp_type_abund.iloc[in_level_2_b[l]][:-1])
								s_type2 = color_picker(samp_type_abund.iloc[in_level_2_b[k]][:-1])
								edge_samp = matchyn(s_type1,s_type2)
								#include the "reverse" edge as well so that cytoscape doesn't have gaps in 
								#it's classification of nodes.
								edge =  [abundance_array['TAXA'][in_level_2_b[l]], abundance_array['TAXA'][in_level_2_b[k]],
															str(adj_matrix_bins[l,k]), str(tot_seen_2_b[l]),str(tot_seen_2_b[k]),
															s_type1[0],s_type2[0],edge_samp[0],s_type1[1],s_type2[1],edge_samp[1],N]
								rev_edge = [abundance_array['TAXA'][in_level_2_b[k]], abundance_array['TAXA'][in_level_2_b[l]],
															str(adj_matrix_bins[l,k]),str(tot_seen_2_b[k]),str(tot_seen_2_b[l]),
															s_type2[0],s_type1[0],edge_samp[0],s_type2[1],s_type1[1],edge_samp[1],N]
								source_target_data_bins += [edge]
								source_target_data_bins += [rev_edge]
	
			
				#turn list into a dataframe
				source_target_frames_bins = pd.DataFrame(source_target_data_bins, columns = ['source','target','weight','source_freq',
									'target_freq','first_sample','second_sample','edge_sample','fscolor','sscolor','edcolor','num_samps'])
				if save:
					flname1 = net_name+'/bins/'+i + '_' + stype+'_adj.tsv'
					flname2 = net_name+'/bins/'+i + '_' + stype+'_list.tsv'
					adjacency_frames_bins.to_csv(flname1, sep = '\t')
					source_target_frames_bins.to_csv(flname2, sep = '\t')
		
		if new_way:
		#Just the matrix product of the normalized abundance vectors will give cosine of the angle between them
		#t1 = time.time()
			#but pearson's correlation coefficient might not make a ton of sense with
			#different sample types. so. let's have the option to average them across sample types.
			pears = array(ab_np_array)
			avgit = False
			if avgit:
				if stype == 'all':
					dsamp_types = diff_samps_types
				else:
					dsamp_types = [stype]
				for idx in dsamp_types:
					selec = lvl_abundance_array.columns.map(lambda x: bool(re.search(idx,x)))
					selllec = where(selec[1:])
					meanss = mean(lvl_abundance_array[lvl_abundance_array.columns[selllec]],axis = 1).values
					vars = std(lvl_abundance_array[lvl_abundance_array.columns[selllec]],axis = 1).values
					pears[:,selllec[0]] = pears[:,selllec[0]] - outer(meanss,ones(len(selllec[0])))
					vars[where(vars == 0)] = 1
					for idx2 in range(len(pears)):
						pears[idx2,selllec[0]] = pears[idx2,selllec[0]]/vars[idx2]
					adj_matrix_pre_pear = (1/pears.shape[1])*dot(pears,transpose(pears))
			else:
				means = mean(pears,axis = 1)
				vars = std(pears,axis =1)
				vars[where(vars == 0)] = 1
				pears = transpose(transpose(pears) - means)
				pears = dot(diag(1/vars),pears)
				adj_matrix_pre_pear = (1/pears.shape[1])*dot(pears,transpose(pears)) - eye(pears.shape[0])
			cutoff = 0.8
			adj_matrix_pre_pear[where(adj_matrix_pre_pear < cutoff)] = 0 
			#print(time.time()-t1)
			degs_p = sum(adj_matrix_pre_pear,1)
			unconnected_p = where(degs_p == 0)
			adj_matrix_pre_pear = delete(adj_matrix_pre_pear,unconnected_p,0)
			adj_matrix_pre_pear = delete(adj_matrix_pre_pear,unconnected_p,1)
			ab_np_array_pear = delete(ab_np_array,unconnected_p,0)
			if ab_np_array_pear.shape[0] != 0:
				in_level_p = delete(in_level,unconnected_p)

				stringency = 0.05
				N = len(ab_np_array[0])
				tot_seen_p = [sum([abd != 0 for abd in row]) for row in ab_np_array]

				ab_masked = mask.masked_values(ab_np_array_pear,0,copy = False)
				the_ns_v = sum(ab_np_array_pear,1)/(ab_masked.min(axis = 1))
				the_ns_v = the_ns_v.data
				the_ns_v = around(the_ns_v)
				the_ns = outer(the_ns_v,ones(len(ab_np_array_pear[0])))
				the_ps_v = sum(ab_np_array_pear,0)/sum(ab_np_array_pear)
				the_ps = outer(ones(len(ab_np_array_pear)),the_ps_v)
				adj_size_p = adj_matrix_pre_pear.shape
				adj_matrix_pear = zeros(adj_size_p) 
				mc_test = mc_pearson(the_ns,the_ps,adj_matrix_pre_pear)#,dsamp_types, lvl_abundance_array.columns)		
				adj_matrix_pear[where(mc_test <stringency)] = adj_matrix_pre_pear[where(mc_test <stringency)]

			
		
				degs2_p = sum(adj_matrix_pear,1)
				unconnected2_p = where(degs2_p == 0)
				adj_matrix_pear = delete(adj_matrix_pear,unconnected2_p,0)
				adj_matrix_pear = delete(adj_matrix_pear,unconnected2_p,1)
				tot_seen_2_p = delete(tot_seen_p,unconnected2_p)
				
				in_level_2_p = delete(in_level_p,unconnected2_p)
			
				adjacency_frames_pear = pd.DataFrame(adj_matrix_pear, index = abundance_array['TAXA'][in_level_2_p], columns = abundance_array['TAXA'][in_level_2_p])
				#create a list of edges - source in one column, target in the other. This is an alternative to the adjacency matrix
				#that might make it easier to load in to cytoscape.
				num_edge_p = count_nonzero(adj_matrix_pear)
				source_target_data_pear = []#empty([num_edge,3], dtype = 'string')	
				for l in range(len(adj_matrix_pear)):
					for k in range(l+1):
							if adj_matrix_pear[l,k] != 0:
								s_type1 = color_picker(samp_type_abund.iloc[in_level_2_p[l]][:-1])
								s_type2 = color_picker(samp_type_abund.iloc[in_level_2_p[k]][:-1])
								edge_samp = matchyn(s_type1,s_type2)
								#include the "reverse" edge as well so that cytoscape doesn't have gaps in 
								#it's classification of nodes.
								edge =  [abundance_array['TAXA'][in_level_2_p[l]], abundance_array['TAXA'][in_level_2_p[k]],
															str(adj_matrix_pear[l,k]), str(tot_seen_2_p[l]),str(tot_seen_2_p[k]),
															s_type1[0],s_type2[0],edge_samp[0],s_type1[1],s_type2[1],edge_samp[1],N]
								rev_edge = [abundance_array['TAXA'][in_level_2_p[k]], abundance_array['TAXA'][in_level_2_p[l]],
															str(adj_matrix_pear[l,k]),str(tot_seen_2_p[k]),str(tot_seen_2_p[l]),
															s_type2[0],s_type1[0],edge_samp[0],s_type2[1],s_type1[1],edge_samp[1],N]
								source_target_data_pear += [edge]
								source_target_data_pear += [rev_edge]
	
		
				#turn list into a dataframe
				source_target_frames_pear = pd.DataFrame(source_target_data_pear, columns = ['source','target','weight','source_freq',
								'target_freq','first_sample','second_sample','edge_sample','fscolor','sscolor','edcolor','num_samps'])
	
				if save:
					flname3 = net_name+'/pears/'+i + '_' + stype+'_adj.tsv'
					flname4 = net_name+'/pears/'+i + '_' + stype+'_list.tsv'
					adjacency_frames_pear.to_csv(flname3, sep = '\t')
					source_target_frames_pear.to_csv(flname4, sep = '\t')
				

#Save adjacency_frames to save the (weighted) adjacency matrix. This can be loaded into 
#cytoscape. Save source_target_frames to save a list of edges. This can also be loaded into
#cytoscape and provides a way to include edge or node attributes.

# print(source_target_frames_pear.keys())
# 
# save = False
# if save:
# 	for j in source_target_frames_pear.keys():
# 		flname1 = net_name+'/bins/'+j+'_adj.tsv'
# 		flname2 = net_name+'/bins/'+j+'_list.tsv'
# 		flname3 = net_name+'/pears/'+j+'_adj.tsv'
# 		flname4 = net_name+'/pears/'+j+'_list.tsv'
# 		adjacency_frames_bins[j].to_csv(flname1, sep = '\t')
# 		source_target_frames_bins[j].to_csv(flname2, sep = '\t')
# 		adjacency_frames_pear[j].to_csv(flname3, sep = '\t')
# 		source_target_frames_pear[j].to_csv(flname4, sep = '\t')


#### TO DO #####################
#
# - Decide how to classify nodes that appear in multiple types of sample at 
#				high frequency
#
#















