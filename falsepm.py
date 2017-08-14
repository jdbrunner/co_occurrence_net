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
	matrices_cor = [contents[j] for j in where([i[-11:] == 'cor_adj.tsv' for i in contents])[0]]
	matrices_thr = [contents[j] for j in where([i[-11:] == 'thr_adj.tsv' for i in contents])[0]]
	#node_tables = [contents[j] for j in where([i[-8:] == 'data.tsv' for i in contents])[0]]
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
	
		
	fposo = True
	fnegs = False
	
	
	networks_cor = dict()
	for mat in matrices_cor:
		undscr = mat.find('_')
		network_adj = pd.read_csv(folder + '/' + mat, sep = '\t')
		networks_cor[mat[undscr+1:-8]] = network_adj
	
	networks_thr = dict()
	for mat in matrices_thr:
		undscr = mat.find('_')
		network_adj = pd.read_csv(folder + '/' + mat, sep = '\t')
		networks_thr[mat[undscr+1:-8]] = network_adj
	
	sort_samples = dict()
	for ky in networks_cor:
		if 'Full' in ky:
			which_samps = where([False,False] + [type(sn) == str for sn in hld_columns_df.columns[2:]])[0]
		else:
			which_samps = where([False,False] + [styp in ky for styp in hld_stypes])[0]
		sort_samples[ky[:-4]] = which_samps
		
# 	for i in sort_samples:
# 		print(i,': ',len(sort_samples[i]))	
# 	sys.exit()
	
	#for each network in the folder, test identification of false +/- on modified holdouts
	#from the network (so only use samples from the correct type). Also test how well the 
	#samples "fit" compared to random data and data from the wrong network
	
	score_results_cor = pd.DataFrame(columns = ['Sample Type'] + [net_type for net_type in networks_cor.keys()])
	score_results_thr = pd.DataFrame(columns = ['Sample Type'] + [net_type for net_type in networks_thr.keys()])


	num_rks = 30
	rank_vs_rank_cor = pd.DataFrame(columns = range(num_rks))
	rank_vs_rank_thr = pd.DataFrame(columns = range(num_rks))
	
	for smp_type in sort_samples:
		print(smp_type)
		cter = 1
		for smp_num in sort_samples[smp_type]:
			print('\r',cter, '/',len(sort_samples[smp_type]), end = '     ')
			cter = cter + 1
			
			the_samp = hld_columns_df.iloc[:,smp_num]
			##We first need to make sure the sample taxa match the network taxa. The sample 
			#could have stuff not in the network
			the_samp.index = hld_columns_df['TAXA']
			nzers = nonzero(the_samp)[0]#find zeros and nonzeros
			zers = where(the_samp == 0)[0] #1D so just take index as int not tuple

			
			if len(nzers)> 1 and sum(the_samp) > 0.0001:
				###Fit Test
				#run the ranking on each network type(correlation)
				fit_scores_cor = pd.DataFrame([[smp_type] + list(zeros(len(networks_cor)))],columns = ['Sample Type'] + [net_type for net_type in networks_cor.keys()])
				for net_type in networks_cor:
					if len(networks_cor[net_type])>1:
						the_sampt = the_samp.loc[networks_cor[net_type]['TAXA']]#now its got only stuff in the network and the right order
						score = ivp_score(networks_cor[net_type], the_sampt)
						fit_scores_cor.loc[0,net_type] = score
				score_results_cor = pd.concat([score_results_cor, fit_scores_cor])
				
				fit_scores_thr = pd.DataFrame([[smp_type] + list(zeros(len(networks_thr)))],columns = ['Sample Type'] + [net_type for net_type in networks_thr.keys()])
				for net_type in networks_thr:
					if len(networks_thr[net_type])>1:
						the_sampt = the_samp.loc[networks_thr[net_type]['TAXA']]#now its got only stuff in the network and the right order
						score = ivp_score(networks_thr[net_type], the_sampt)
						fit_scores_thr.loc[0,net_type] = score
				score_results_thr = pd.concat([score_results_thr, fit_scores_thr])
			
				####### Trim to full network and sort
				the_samp_cor_trimmed = the_samp.loc[networks_cor['Full_cor']['TAXA']]#now its got only stuff in the network and the right order
				the_samp_thr_trimmed = the_samp.loc[networks_thr['Full_thr']['TAXA']]#now its got only stuff in the network and the right order
				
				smp_srted_cor = the_samp_cor_trimmed.sort_values(ascending=False)
				smp_srted_thr = the_samp_thr_trimmed.sort_values(ascending=False)
				
				rking_results_cor = pd.DataFrame([zeros(num_rks)], columns = range(num_rks))
				rking_results_thr = pd.DataFrame([zeros(num_rks)], columns = range(num_rks))
				
				####### Remove the nth ranked by abundance and see where it comes in the ranking from IVP
				for tax_num in range(num_rks):
					ind_loc_cor = where(the_samp_cor_trimmed.index == smp_srted_cor.index[tax_num])
					ind_loc_thr = where(the_samp_thr_trimmed.index == smp_srted_thr.index[tax_num])
					fneg_samp_cor = array(the_samp_cor_trimmed.values)
					fneg_samp_cor[ind_loc_cor] = 0
					fneg_samp_thr = array(the_samp_thr_trimmed.values)
					fneg_samp_thr[ind_loc_thr] = 0
					this_rking_cor = diffusion_ivp([],[],networks_cor['Full_cor'].values[:,1:] , sample = fneg_samp_cor, all = True)[0]
					this_rking_thr = diffusion_ivp([],[],networks_thr['Full_thr'].values[:,1:] , sample = fneg_samp_thr, all = True)[0]
					the_rk_cor = where(array(flat_two_deep(this_rking_cor)) == ind_loc_cor[0][0])[0][0]
					the_rk_thr = where(array(flat_two_deep(this_rking_thr)) == ind_loc_thr[0][0])[0][0]
					rking_results_cor.loc[0,tax_num] = the_rk_cor - len(nzers) + 1
					rking_results_thr.loc[0,tax_num] = the_rk_thr - len(nzers) + 1
				
				rank_vs_rank_cor = pd.concat([rank_vs_rank_cor, rking_results_cor])
				rank_vs_rank_thr = pd.concat([rank_vs_rank_thr, rking_results_thr])
		
		print('\n')		
				
				#### Here we could add permutation based null model
	
	score_results_cor.index = range(len(score_results_cor))
	score_results_thr.index = range(len(score_results_thr))
	
	### Fit scores for randomly generated samples
	print('random samples')
	random_samples = make_sample(hld_columns_df, numb = len(hld_columns_df.columns) - 2)
	
		
	rand_score_results_cor = pd.DataFrame(columns = ['Sample Type'] + [net_type for net_type in networks_cor.keys()])
	rand_score_results_thr = pd.DataFrame(columns = ['Sample Type'] + [net_type for net_type in networks_thr.keys()])
	
	cter = 1
	for smp_num in range(2,len(random_samples.columns)):	
		print('\r',cter, '/',len(random_samples.columns) - 2, end = '     ')
		cter = cter + 1		
		the_samp = random_samples.iloc[:,smp_num]
		##We first need to make sure the sample taxa match the network taxa. The sample 
		#could have stuff not in the network
		the_samp.index = random_samples['TAXA']
		nzers = nonzero(the_samp)[0]#find zeros and nonzeros
		zers = where(the_samp == 0)[0] #1D so just take index as int not tuple
		
		if len(nzers)>1:
			###Fit Test
			#run the ranking on each network type(correlation)
			fit_scores_cor = pd.DataFrame([['random'] + list(zeros(len(networks_cor)))],columns = ['Sample Type'] + [net_type for net_type in networks_cor.keys()])
			for net_type in networks_cor:
				if len(networks_cor[net_type])>1:
					the_sampt = the_samp.loc[networks_cor[net_type]['TAXA']]#now its got only stuff in the network and the right order
					score = ivp_score(networks_cor[net_type], the_sampt)
					fit_scores_cor.loc[0,net_type] = score
			rand_score_results_cor = pd.concat([rand_score_results_cor, fit_scores_cor])
			
			fit_scores_thr = pd.DataFrame([['random'] + list(zeros(len(networks_thr)))],columns = ['Sample Type'] + [net_type for net_type in networks_thr.keys()])
			for net_type in networks_thr:
				if len(networks_thr[net_type])>1:
					the_sampt = the_samp.loc[networks_thr[net_type]['TAXA']]#now its got only stuff in the network and the right order
					score = ivp_score(networks_thr[net_type], the_sampt)
					fit_scores_thr.loc[0,net_type] = score
			rand_score_results_thr = pd.concat([rand_score_results_thr, fit_scores_thr])
		
	
	
	############## Plotting
	print('\n plotting')
	
	plt_folder = folder + '/validation_plots'
	
	#### Histograms of fit scores
	fit_histos_cor, axarr_cor = plt.subplots(2, 1, figsize = (15,10))
	axarr_cor[0].hist(array(score_results_cor.Full_cor.values), bins = 50)
	axarr_cor[0].set_title('Held Out Columns', fontsize=10)
	axarr_cor[1].hist(array(rand_score_results_cor.Full_cor.values), bins = 50)
	axarr_cor[1].set_title('Random Columns', fontsize=10)
	fit_histos_cor.suptitle('Full Correlation Network Fit', fontsize = 20)
	
	fit_histos_cor.savefig(plt_folder + '/fit_histo_cor.png')
	plt.close(fit_histos_cor)
	
	fit_histos_thr, axarr_thr = plt.subplots(2, 1, figsize = (15,10))
	axarr_thr[0].hist(array(score_results_thr.Full_thr.values), bins = 50)
	axarr_thr[0].set_title('Held Out Columns', fontsize=10)
	axarr_thr[1].hist(array(rand_score_results_thr.Full_thr.values), bins = 50)
	axarr_thr[1].set_title('Random Columns', fontsize=10)
	fit_histos_thr.suptitle('Full Threshhold Network Fit', fontsize = 20)
	
	fit_histos_thr.savefig(plt_folder + '/fit_histo_thr.png')
	plt.close(fit_histos_thr)

	
	###### Bar chart of correct identification
	score_res_nofull_cor = score_results_cor.drop('Full_cor', axis = 1)
	score_res_nofull_thr = score_results_thr.drop('Full_thr', axis = 1)
	max_cor = score_res_nofull_cor.iloc[:,1:].idxmax(axis = 1)
	max_thr = score_res_nofull_thr.iloc[:,1:].idxmax(axis = 1)
	score_res_nofull_cor['Winner'] = max_cor
	score_res_nofull_thr['Winner'] = max_thr
	
	
	perc_right_cor = []
	for typ in unique(score_res_nofull_cor['Sample Type']):
		jj =score_res_nofull_cor[score_res_nofull_cor['Sample Type'] == typ].index[0]
		correct = sum([int(typ in score_res_nofull_cor.loc[j,'Winner']) for j in score_res_nofull_cor[score_res_nofull_cor['Sample Type'] == typ].index])
		tots = len(score_res_nofull_cor[score_res_nofull_cor['Sample Type'] == typ])
		perc_right_cor =perc_right_cor + [correct/tots]

	perc_right_thr = []
	for typ in unique(score_res_nofull_thr['Sample Type']):
		correct = sum([int(typ in score_res_nofull_thr.loc[j,'Winner']) for j in score_res_nofull_thr[score_res_nofull_thr['Sample Type'] == typ].index])
		tots = len(score_res_nofull_thr[score_res_nofull_thr['Sample Type'] == typ])
		perc_right_thr =perc_right_thr + [correct/tots]
		
	score_res_nofull_cor.drop('Winner', axis = 1, inplace = True)
	score_res_nofull_thr.drop('Winner', axis = 1, inplace = True)
	
	win_perc, ax_win = plt.subplots(2,1,tight_layout=True, figsize = (15,10))
	ax_win[0].bar(range(len(perc_right_cor)), perc_right_cor)
	ax_win[0].set_title('Correlation Network', fontsize = 10)
	ax_win[0].set_xticks(range(len(unique(score_res_nofull_cor['Sample Type']))))
	ax_win[0].set_xticklabels(unique(score_res_nofull_cor['Sample Type']))
	labels = ax_win[0].get_xticklabels()
	plt.setp(labels, rotation=40, fontsize=8)
	ax_win[1].bar(range(len(perc_right_thr)), perc_right_thr)
	ax_win[1].set_title('Co-occurrence (Threshhold) Network', fontsize = 10)
	ax_win[1].set_xticks(range(len(unique(score_res_nofull_thr['Sample Type']))))
	ax_win[1].set_xticklabels(unique(score_res_nofull_thr['Sample Type']))
	labels = ax_win[1].get_xticklabels()
	plt.setp(labels, rotation=40, fontsize=8)	
	win_perc.savefig(plt_folder + '/correct_ids.png')
	plt.close(win_perc)
	
	
	##### Bar chart of average fit for sample type by network
	avg_scores_cor = dict()
	avg_scores_thr = dict()
	for typ in unique(score_res_nofull_cor['Sample Type']):
		avg_scores_cor[typ] = mean(score_res_nofull_cor[score_res_nofull_cor['Sample Type'] == typ].iloc[1:], axis = 0)
	
	for typ in unique(score_res_nofull_thr['Sample Type']):
		avg_scores_thr[typ] = mean(score_res_nofull_thr[score_res_nofull_thr['Sample Type'] == typ].iloc[1:], axis = 0)
	
	
	for typ in unique(score_res_nofull_cor['Sample Type']):
		av, ax = plt.subplots(1,1, tight_layout=True)
		ax.bar(range(len(avg_scores_cor[typ])), avg_scores_cor[typ])
		ax.set_title(typ + ' Fit Scores, Correlation')
		ax.set_xticks(range(len(score_res_nofull_cor.columns[1:])))
		ax.set_xticklabels(score_res_nofull_cor.columns[1:])
		labels = ax.get_xticklabels()
		plt.setp(labels, rotation=60, fontsize=10)
		av.savefig(plt_folder + '/fit_score_cor_'+typ+'.png')
		plt.close(av)

	for typ in unique(score_res_nofull_thr['Sample Type']):
		av, ax = plt.subplots(1,1, tight_layout=True)
		ax.bar(range(len(avg_scores_thr[typ])), avg_scores_thr[typ])
		ax.set_title(typ + ' Fit Scores, Threshhold')
		ax.set_xticks(range(len(score_res_nofull_thr.columns[1:])))
		ax.set_xticklabels(score_res_nofull_thr.columns)
		labels = ax.get_xticklabels()
		plt.setp(labels, rotation=60, fontsize=10)
		av.savefig(plt_folder + '/fit_score_thr_'+typ+'.png')
		plt.close(av)
		
	
	##### Scatter plot of Diffusion Rk vs Abundance Rk	
	
	scattr_cor = []
	rank_vs_rank_cor.index = range(len(rank_vs_rank_cor))
	for i in rank_vs_rank_cor.columns:
		scattr_cor = scattr_cor + [(i, rank_vs_rank_cor.loc[j,i]) for j in rank_vs_rank_cor.index]
	
	scattr_thr = []
	rank_vs_rank_thr.index = range(len(rank_vs_rank_thr))
	for i in rank_vs_rank_thr.columns:
		scattr_thr = scattr_thr + [(i, rank_vs_rank_thr.loc[j,i]) for j in rank_vs_rank_thr.index]
	
# 	print(scattr_cor)
	sc_cor, ax_cor = plt.subplots(1,1)
	ax_cor.scatter(*zip(*scattr_cor))
	ax_cor.set_title('Rank by Abundance vs Rank by Diffusion, Correlation')
	ax_cor.set_xlabel('Rank by Abundance')
	ax_cor.set_ylabel('Rank by Diffusion')
	
	sc_cor.savefig(plt_folder + '/rank_v_rank_cor.png')
	plt.close(sc_cor)
	
	sc_thr, ax_thr = plt.subplots(1,1)
	ax_thr.scatter(*zip(*scattr_thr))
	ax_thr.set_title('Rank by Abundance vs Rank by Diffusion, Correlation')
	ax_thr.set_xlabel('Rank by Abundance')
	ax_thr.set_ylabel('Rank by Diffusion')
	
	sc_thr.savefig(plt_folder + '/rank_v_rank_thr.png')
	plt.close(sc_thr)
	
	
	
	
	
	
	