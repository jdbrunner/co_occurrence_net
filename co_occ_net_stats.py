#########################################################################################
#																						#		
#				Co-Occurrence Network Stats												#
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

#import networkx as nx

#####Variables
##User input, should be entered as an argument
folder = sys.argv[1]
flder = os.listdir(folder)
matrices = [flder[j] for j in where([i[-7:] == 'adj.tsv' for i in flder])[0]]
node_table = [flder[j] for j in where([i[-8:] == 'data.tsv' for i in flder])[0]]
numnets = len(matrices)







		
######Import the network as list of edges
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 


for that_net in range(numnets):
	node_data = pd.read_csv(folder + '/' + node_table[that_net], sep = '\t')
	node_data.index = node_data['TAXA']
	network_mat = pd.read_csv(folder + '/' + matrices[that_net],sep = "\t")
	if len(node_data) < 10:
		continue
		



	adj_mat = network_mat.values[:,1:]
	communities = com_clust(network_mat)
	spect = spectral_cluster(adj_mat)
	spectcols = color_picker2(unique(spect))
	commcols = color_picker2(unique(communities[0]))
	
	node_data['spect_cluster'] = zeros(len(node_data))
	node_data['commun_cluster'] = zeros(len(node_data))
	node_data['spect_color'] = zeros(len(node_data))
	node_data['comm_color'] = zeros(len(node_data))
	spect_colors = [spectcols[i] for i in spect]
	comm_colors = [commcols[i] for i in communities[0]]
	for j in range(len(node_data)):
		node_data.loc[network_mat.iloc[j]['TAXA'],'spect_cluster'] = spect[j]
		node_data.loc[network_mat.iloc[j]['TAXA'],'commun_cluster'] = communities[0][j]
		node_data.loc[network_mat.iloc[j]['TAXA'],'spect_color'] = spect_colors[j]
		node_data.loc[network_mat.iloc[j]['TAXA'],'comm_color'] = comm_colors[j]
	
	node_data.drop('Unnamed: 0', axis = 1, inplace = True)
	node_data.drop('TAXA', axis = 1, inplace = True)
	node_data.to_csv(folder + '/' + node_table[that_net], sep = '\t')
















