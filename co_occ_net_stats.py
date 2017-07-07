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
edges = [flder[j] for j in where([i[-8:] == 'list.tsv' for i in flder])[0]]
numnets = len(matrices)







		
######Import the network as list of edges
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 


for that_net in range(numnets):
	network_edges = pd.read_csv(folder + '/' + edges[that_net], sep = '\t')
	network_edges.index = network_edges['source']
	network_mat = pd.read_csv(folder + '/' + matrices[that_net],sep = "\t")
	if len(network_edges) < 10:
		continue
		



	adj_mat = network_mat.values[:,1:]
	communities = com_clust(network_mat)
	spect = spectral_cluster(adj_mat)
	spectcols = color_picker2(unique(spect))
	commcols = color_picker2(unique(communities[0]))
	network_edges['spect_cluster'] = zeros(len(network_edges))
	network_edges['commun_cluster'] = zeros(len(network_edges))
	network_edges['spect_color'] = zeros(len(network_edges))
	network_edges['comm_color'] = zeros(len(network_edges))
	spect_colors = [spectcols[i] for i in spect]
	comm_colors = [commcols[i] for i in communities[0]]
	for j in range(len(network_mat)):
		network_edges.loc[network_mat.iloc[j]['TAXA'],'spect_cluster'] = spect[j]
		network_edges.loc[network_mat.iloc[j]['TAXA'],'commun_cluster'] = communities[0][j]
		network_edges.loc[network_mat.iloc[j]['TAXA'],'spect_color'] = spect_colors[j]
		network_edges.loc[network_mat.iloc[j]['TAXA'],'comm_color'] = comm_colors[j]

	end_name = edges[that_net].find('.tsv')
	cled_name = 'clustered/' + edges[that_net][:end_name] + '_clustered.tsv'
	network_edges.to_csv(folder + '/' + cled_name, sep = '\t')
















