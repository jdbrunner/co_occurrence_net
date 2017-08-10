#########################################################################################
#																						#		
#				False positives & negatives, network fit								#
#																						#
#										Jim Brunner 									#
#																						#
#																						#
#																						#
#########################################################################################
#
# We want test the idea that diffusion can be used to identify possible false positives
# or false negatives
#
# To do this, we can use hold out columns of the data that we modify by removing or adding 
#taxa, and check whether the method identifies the changes.
#
# We also want to see that samples of hold out data fit a network better than random samples
# so we do that too
#
#
## modules needed
from pylab import *
import pandas as pd
###pandas DataFrame usage: df['column_name'] or df.loc['row_name' (index),'column'] or df.iloc[row_number, 'column']
import sys
import matplotlib.pyplot as plt
from co_occ_funs import *
import os
from random import sample

if __name__ == '__main__':
	abund_name = sys.argv[1]
	folder = sys.argv[2]
	contents = os.listdir(folder)
	#list the adjacency matrices, the node data, and the holdouts.
	matrices = [contents[j] for j in where([i[-7:] == 'adj.tsv' for i in contents])[0]]
	node_tables = [contents[j] for j in where([i[-8:] == 'data.tsv' for i in contents])[0]]
	held_outs = [contents[j] for j in where([i[-8:] == 'held.tsv' for i in contents])[0]][0]
	
	held_out_df = pd.read_csv(folder + '/' + held_outs, index_col = 0, sep = '\t')
	abund_df = pd.read_csv(abund_name, sep = ' ')
	hld_inds = [0,1] + list(held_out_df['Holdouts'])
	hld_columns_df = abund_df.iloc[:,hld_inds]
	lvs = unique(hld_columns_df['LEVEL'].values)
	# the level name is the name of a folder in the path
	this_lvl = lvs[where(['/'+nm+'/' in folder for nm in lvs])][0]
	#now we can trim off all the data from the wrong taxonomic level
	hld_columns_df = hld_columns_df.loc[where(hld_columns_df['LEVEL'] == this_lvl)]
	hld_stypes = array([smp[:-10] for smp in hld_columns_df.columns[2:]])
	for sty in range(len(hld_stypes)):
		if 'female' in hld_stypes[sty]:
			hld_stypes[sty] =hld_stypes[sty][:-7]
		elif 'male' in hld_stypes[sty]:
			hld_stypes[sty] = hld_stypes[sty][:-5]
	
	all_fit_scores_samp = []
	all_fit_scores_rand = []
	all_fit_scores_wrong = []
	
	all_fn_cc_rks = []
	all_fn_ov_rks = []
	all_fp_cc_rks = []
	all_fp_ov_rks = []
	
	
	con = 0.2
	
	fposo = True
	fnegs = False
	
	#for each network in the folder, test identification of false +/- on modified holdouts
	#from the network (so only use samples from the correct type). Also test how well the 
	#samples "fit" compared to random data and data from the wrong network
	for mat in matrices:
		#For each network, I'll make the following lists:
		#-fit of hold out columns
		#-fit of random samples
		#-fit of samples from the wrong type (not for full)
		#And the following lists of lists:
		#-false negative connected component ranks (a list for each holdout)
		#-false negative overall ranks (a list for each holdout)
		#-false postive connected component ranks (a list for each holdout)
		#-false postive overall ranks (a list for each holdout)
		print(mat)
		network_adj = pd.read_csv(folder + '/' + mat, sep = '\t')
		if len(network_adj.values)>1:
			#identify which heldout samples go with the network in question
			if 'Full' in mat:
				which_samps = where([False,False] + [type(sn) == str for sn in hld_columns_df.columns[2:]])[0]
				not_these = array([])
			else:
				which_samps = where([False,False] + [styp in mat for styp in hld_stypes])[0]
				not_these = where([False,False] + [not(styp in mat) for styp in hld_stypes])[0]
		
			
			fit_scores_samps = []
		
			ccrks_false_negs = []
			ccrks_false_pos = []
			overrks_false_negs = []
			overrks_false_pos = []
		
			for hld in which_samps:
				the_samp = hld_columns_df.iloc[:,hld]
				##We first need to make sure the sample taxa match the network taxa. The sample 
				#could have stuff not in the network
				the_samp.index = hld_columns_df['TAXA']
				the_samp = the_samp.loc[network_adj['TAXA']]#now its got only stuff in the network and the right order
				nzers = nonzero(the_samp)[0]#find zeros and nonzeros
				zers = where(the_samp == 0)[0] #1D so just take index as int not tuple
				###Fit Test
				#we run the ranking procedure
				samp_tot = sum(the_samp.values)
				if len(nzers)>1:
					ranking1 = diffusion_ivp([],[], network_adj.values[:,1:], sample = the_samp.values, all = True)[0]
					samp_orded = the_samp.values[flat_two_deep(ranking1)]
					geoms = con**array(range(len(samp_orded)))
					score = dot(geoms, samp_orded/samp_tot)
# 					if score == 1:
# 						print(samp_tot)
# 						print(samp_orded)
					fit_scores_samps = fit_scores_samps + [score]
			
			
				#### False +/-
				#make a modified holdout
				if samp_tot !=0:
					if fnegs:
						nfnegs = int(ceil(0.1*len(nzers)))#choose how many FN and FP
					else:
						nfnegs = 0
					if fposo:
						nfpos = int(ceil(0.1*len(zers)))
					else:
						nfpos = 0
					false_negs = nzers[sample(range(len(nzers)), k = nfnegs)]#choose FNs and FPs
					false_pos = zers[sample(range(len(zers)),k = nfpos)]
					the_samp.iloc[false_negs] = zeros(len(false_negs))
					#make the false positive a uniform random times the average nonzero
					nzavg = mean(the_samp[nzers])
					the_samp.iloc[false_pos] = rand(len(false_pos))*nzavg
					#now run the ranking procedure. The result is a list of the indices of nodes in the graph
					#(thus taxa). 
					ranking2 = diffusion_ivp([],[], network_adj.values[:,1:], sample = the_samp.values, all = True)[0]
					#put everything into an inner list
					ranking2 = [[num] if not(isinstance(num,list)) else num for num in ranking2]
					#number connected components?
					num_ccs = len(ranking2)
					#flatten ranking
					fltrank2 = flat_two_deep(ranking2)
					#flatten each cc so its more easily searched
					cc_flt_rk = [flat_one_deep(cco) for cco in ranking2]
					###where did the false negs end up ranked?
					fn_ccrks = []
					fn_ovrrks = []
					fp_ccrks = []
					fp_ovrrks = []
					for fn in false_negs:
						#rank of connected component
						ccrk = argwhere([fn in cc for cc in cc_flt_rk])
						comrk = ccrk[0]
						#and overall rank
						overrk = argwhere([fn == nd for nd in fltrank2])
						#save them
						fn_ccrks = fn_ccrks + [comrk]
						fn_ovrrks = fn_ovrrks + [overrk]
					for fp in false_pos:
						#rank of connected component
						ccrk = argwhere([fp in cc for cc in cc_flt_rk])
						comrk = ccrk[0]
						#and overall rank
						overrk = argwhere([fp == nd for nd in fltrank2])
						#save them
						fp_ccrks = fp_ccrks + [comrk]
						fp_ovrrks = fp_ovrrks + [overrk]
			

					ccrks_false_negs = ccrks_false_negs + [rk/num_ccs for rk in fn_ccrks]
					ccrks_false_pos = ccrks_false_pos + [rk/num_ccs for rk in fp_ccrks]
					overrks_false_negs = overrks_false_negs + [rk/len(the_samp) for rk in fn_ovrrks]
					overrks_false_pos = overrks_false_pos + [rk/len(the_samp) for rk in fp_ovrrks]
# 		

			fit_scores_rands = []
			for rs in range(100):
				rsamp = make_sample(hld_columns_df)			 
				##We first need to make sure the sample taxa match the network taxa. The sample 
				#could have stuff not in the network
				rsamp.index = hld_columns_df['TAXA']
				rsamp = rsamp.loc[network_adj['TAXA']]#now its got only stuff in the network and the right order
				nzers = nonzero(rsamp['abundance'])[0]#find zeros and nonzeros
				###Fit Test
				#we run the ranking procedure
				samp_tot = sum(rsamp['abundance'].values)
				if len(nzers) > 1:
					ranking1 = diffusion_ivp([],[], network_adj.values[:,1:], sample = rsamp['abundance'].values, all = True)[0]
					samp_orded = rsamp['abundance'].values[flat_two_deep(ranking1)]
					geoms = con**array(range(len(samp_orded)))
					score = dot(geoms, samp_orded/samp_tot)
					fit_scores_rands = fit_scores_rands + [score]
					
			
			fit_scores_wrong = []
			for hld in not_these:
				the_samp = hld_columns_df.iloc[:,hld]
				##We first need to make sure the sample taxa match the network taxa. The sample 
				#could have stuff not in the network
				the_samp.index = hld_columns_df['TAXA']
				the_samp = the_samp.loc[network_adj['TAXA']]#now its got only stuff in the network and the right order
				nzers = nonzero(the_samp)[0]#find zeros and nonzeros
				###Fit Test
				#we run the ranking procedure
				samp_tot = sum(the_samp.values)
				if len(nzers)>1:
					ranking1 = diffusion_ivp([],[], network_adj.values[:,1:], sample = the_samp.values, all = True)[0]
					samp_orded = the_samp.values[flat_two_deep(ranking1)]
					geoms = con**array(range(len(samp_orded)))
					score = dot(geoms, samp_orded/samp_tot)
# 					if score == 1:
# 						print(samp_tot)
# 						print(samp_orded)
					fit_scores_wrong = fit_scores_wrong + [score]		
	
	
			all_fit_scores_samp = all_fit_scores_samp + fit_scores_samps
			all_fit_scores_rand = all_fit_scores_rand + fit_scores_rands
			all_fit_scores_wrong = all_fit_scores_wrong + fit_scores_wrong
	
			all_fn_cc_rks = all_fn_cc_rks + [y for x in ccrks_false_negs for y in x]
			all_fn_ov_rks = all_fn_ov_rks + [y for x in overrks_false_negs for y in x]
			all_fp_cc_rks = all_fp_cc_rks + [y for x in ccrks_false_pos for y in x]
			all_fp_ov_rks = all_fp_ov_rks + [y for x in overrks_false_pos for y in x]
			
	
	fit_histos, axarr = plt.subplots(3, 1)
	axarr[0].hist(array(all_fit_scores_samp), bins = 50)
	axarr[0].set_title('Held Out Columns', fontsize=10)
	axarr[1].hist(array(all_fit_scores_rand), bins = 50)
	axarr[1].set_title('Random Columns', fontsize=10)
	axarr[2].hist(array(all_fit_scores_wrong), bins = 50)
	axarr[2].set_title('Wrong Network Columns', fontsize=10)
	fit_histos.suptitle('Fit Scores', fontsize = 20)
		
	
	rank_histos, axarr2 = plt.subplots(2,1)
	rank_histos_pos, axarr3 = plt.subplots(2,1)
	
	axarr2[0].hist(array(all_fn_cc_rks), bins = 25)
	axarr2[0].set_title('False Negative Connected Component Ranks', fontsize=10)
	axarr3[0].hist(array(all_fp_cc_rks), bins = 25)
	axarr3[0].set_title('False Positive Connected Component Ranks', fontsize=10)
	axarr2[1].hist(array(all_fn_ov_rks), bins = 25)
	axarr2[1].set_title('False Negative Overall Ranks', fontsize=10)
	axarr3[1].hist(array(all_fp_ov_rks), bins = 25)
	axarr3[1].set_title('False Postive Overall Ranks', fontsize=10)
	rank_histos.suptitle('Rank of Modified Nodes', fontsize=20)
	rank_histos_pos.suptitle('Rank of Modified Nodes', fontsize=20)

	show()
	
	
	
	
	
	