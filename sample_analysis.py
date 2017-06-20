#########################################################################################
#																						#		
#				Sample Analysis															#
#																						#
#										Jim Brunner 									#
#																						#
#																						#
#																						#
#########################################################################################
#
# We want to use the co-occurence network to say something about a
#	sample
#
#
## modules needed
from pylab import *
import pandas as pd
import sys
import matplotlib.pyplot as plt

######################
def est_prob(source,induced,whole,N):
	'''Use distribution factorization from random markov field to estimate the probability of
	seeing the subset of taxa in the sample. Takes in the ROW NUMBERS of the represented edges,
	 along with the whole network.'''
	##use the whole network to get a normalizing constant...on the other hand for many purposes
	#this doesn't matter....
	num_tot_edges = int(len(whole))
	all_freq = array([whole['source_freq'][j] for j in range(0,num_tot_edges)])
	s_freq = delete(all_freq,range(1,num_tot_edges,2))
	t_freq = delete(all_freq,range(0,num_tot_edges,2))
	all_st_freq = array([whole['weight'][j] for j in range(0,num_tot_edges)])
	all_p_s_given_t = all_st_freq/all_freq
	p_s_given_t = delete(all_p_s_given_t,range(1,num_tot_edges,2))
	p_t_given_s = delete(all_p_s_given_t,range(0,num_tot_edges,2))
	normalizer = ones(num_tot_edges//2) + s_freq/N + t_freq/N + 0.5*p_s_given_t + 0.5*p_t_given_s
	Z = prod(normalizer)
	#now use the ROW NUMBERS OF THE represented edges and present edges
	edge_there = argwhere([node in induced[0] for node in source[0]])
	repped = delete(source,edge_there)
	uni_induced = delete(induced,range(1,len(induced[0]),2))
	print(len(uni_induced))
	print(len(repped))
	psi1 = log([0.5*all_p_s_given_t[i] + 0.5*all_p_s_given_t[i+1] for i in uni_induced])
	psi2 = log([all_freq[i]/N for i in repped])
	produ = sum(psi1) + sum(psi2) - log(Z)
	return produ




######################
sample_name = sys.argv[1]
coin_net_name = sys.argv[2]
coin_mat_name = sys.argv[3]
coocc_net_name = sys.argv[4]
coocc_mat_name = sys.argv[5]
level = sys.argv[6]

sample = pd.read_csv(sample_name, sep = ' ')
coin_net = pd.read_csv(coin_net_name, sep = '\t')
coin_mat = pd.read_csv(coin_mat_name, sep = '\t', index_col = 0)
coocc_net = pd.read_csv(coocc_net_name, sep = '\t')
coocc_mat = pd.read_csv(coocc_mat_name, sep = '\t', index_col = 0)

#trim the sample to just the level we are looking at
sample = sample.iloc[where([i == level for i in sample['LEVEL']])]

#who is there
N = len(sample)
present = nonzero(sample['abundance'])[0]
present_names = sample['TAXA'][present]
print(len(present))

#construct the network from the sample.
# samp_net = zeros([N,N])
# for ind in present:
# 	samp_net[ind] = 1
# samp_net = samp_net + transpose(samp_net)
# sample_df = pd.DataFrame(samp_net, columns = sample['TAXA'], index = sample['TAXA'])

#get the induced subgraph of the coincidence and cooccurrence networks.
# induced_coin = coin_mat[present_names]
# induced_coin = induced_coin.loc[present_names]
induced_edges_coin = where([(coin_net['source'][i] in array(present_names)) & (coin_net['target'][i] in array(present_names))
			for i in range(len(coin_net))])
induced_coin_net = coin_net.iloc[induced_edges_coin]

# induced_coocc = coocc_mat[present_names]
# induced_coocc = induced_coocc.loc[present_names]
induced_edges_coocc = where([(coocc_net['source'][i] in array(present_names)) & (coocc_net['target'][i] in array(present_names))
			for i in range(len(coocc_net))])
repped_sources = where([coocc_net['source'][i] in array(present_names) for i in range(len(coocc_net))])
samp_prob = est_prob(repped_sources,induced_edges_coocc,coocc_net,coocc_net['num_samps'][0])
induced_coocc_net = coocc_net.iloc[induced_edges_coocc]
print(samp_prob)

#save the induced networks for loading into cytoscape.
# induced_coin_net.to_csv(sample_name[:-4]+'_induced_coin.tsv', sep = '\t')
# induced_coocc_net.to_csv(sample_name[:-4]+'_induced_coocc.tsv', sep = '\t')


 