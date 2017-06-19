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
induced_coocc_net = coocc_net.iloc[induced_edges_coocc]

induced_coin_net.to_csv(sample_name[:-4]+'_induced_coin.tsv', sep = '\t')
induced_coocc_net.to_csv(sample_name[:-4]+'_induced_coocc.tsv', sep = '\t')

 