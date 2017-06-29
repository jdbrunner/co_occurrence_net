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
coin_matrices = [flder[j] for j in where([i[-12:] == 'coin_adj.tsv' for i in flder])[0]]
coin_edges = [flder[j] for j in where([i[-13:] == 'coin_list.tsv' for i in flder])[0]]
coocc_matrices = [flder[j] for j in where([i[-13:] == 'coocc_adj.tsv' for i in flder])[0]]
coocc_edges = [flder[j] for j in where([i[-14:] == 'coocc_list.tsv' for i in flder])[0]]
numnets_coocc = len(coocc_matrices)
numnets_coin = len(coin_matrices)

# csv_mat_name = sys.argv[1]
# csv_edges_name = sys.argv[2]


#edges = False


	
###### Number of Edges vs Random ############
##
# - Compare number of types of edges vs expected number of that type in a random graph
#   - Use ER random graph
#	- For edges within a sample type this just number of edges vs E(number of edges in 
# 		an ER random graph with number of nodes = number of nodes of that sample type)
#	- For intertype edges this will be harder. I'll have to check this but it should be
#		E(edges of ER random graph) - \sum_{types} E(edges of subgraph)
# 





		
######Import the network as list of edges
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 

for this_net in range(numnets_coin):
	network_edges = pd.read_csv(folder + '/' + coin_edges[this_net], sep = '\t')
	network_edges.index = network_edges['source']
	network_mat = pd.read_csv(folder + '/' + coin_matrices[this_net],sep = "\t")
	if len(network_edges) < 10:
		continue
	#	net_graph = nx.from_pandas_dataframe(network_edges, 'source', 'target', 'weight')

	edge_an = False
	if edge_an:
		these_types = sys.argv[3]
		if ',' in these_types:
			commas = where([letter == ',' for letter in these_types])
			commas = append(commas,len(these_types)-1)
			commas = insert(commas,0,0)
			the_types = []
			for i in range(1,len(commas)):
				the_types = the_types + ["".join(these_types[commas[i-1]+1:commas[i]])]
			these_types = the_types
		else:
			these_types = ["".join(these_types[1:-1])]
		
		num_edges = len(network_edges)

		sample_types = unique(network_edges['edge_sample'])

		mean_weights = [mean(network_edges['weight'].loc[where(network_edges['edge_sample'] == type)])
								for type in sample_types]

		edge_props = [float(len(network_edges['weight'].loc[where(network_edges['edge_sample'] == type)]))/float(num_edges)
								for type in sample_types]
		####### Just getting mean weights of edges within a sample type (and intertype), and 
		####proportion of said edges. And then making a bar plot.
		mean_weights_dict = dict()
		for i in range(len(sample_types)):
			mean_weights_dict[sample_types[i]] = mean_weights[i]

		edge_props_dict = dict()
		for i in range(len(sample_types)):
			edge_props_dict[sample_types[i]] = edge_props[i]
	

		e_prob = edge_prob(network_edges)


	clus_an = True
	if clus_an:
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

	end_name = coin_edges[this_net].find('.tsv')
	cled_name = 'clustered/' + coin_edges[this_net][:end_name] + '_clustered.tsv'
	network_edges.to_csv(folder + '/' + cled_name, sep = '\t')


for that_net in range(numnets_coocc):
	network_edges = pd.read_csv(folder + '/' + coocc_edges[that_net], sep = '\t')
	network_edges.index = network_edges['source']
	network_mat = pd.read_csv(folder + '/' + coocc_matrices[that_net],sep = "\t")
	if len(network_edges) < 10:
		continue
		
	#	net_graph = nx.from_pandas_dataframe(network_edges, 'source', 'target', 'weight')

	edge_an = False
	if edge_an:
		these_types = sys.argv[3]
		if ',' in these_types:
			commas = where([letter == ',' for letter in these_types])
			commas = append(commas,len(these_types)-1)
			commas = insert(commas,0,0)
			the_types = []
			for i in range(1,len(commas)):
				the_types = the_types + ["".join(these_types[commas[i-1]+1:commas[i]])]
			these_types = the_types
		else:
			these_types = ["".join(these_types[1:-1])]
		
		num_edges = len(network_edges)

		sample_types = unique(network_edges['edge_sample'])

		mean_weights = [mean(network_edges['weight'].loc[where(network_edges['edge_sample'] == type)])
								for type in sample_types]

		edge_props = [float(len(network_edges['weight'].loc[where(network_edges['edge_sample'] == type)]))/float(num_edges)
								for type in sample_types]
		####### Just getting mean weights of edges within a sample type (and intertype), and 
		####proportion of said edges. And then making a bar plot.
		mean_weights_dict = dict()
		for i in range(len(sample_types)):
			mean_weights_dict[sample_types[i]] = mean_weights[i]

		edge_props_dict = dict()
		for i in range(len(sample_types)):
			edge_props_dict[sample_types[i]] = edge_props[i]
	

		e_prob = edge_prob(network_edges)


	clus_an = True
	if clus_an:
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

	end_name = coocc_edges[that_net].find('.tsv')
	cled_name = 'clustered/' + coocc_edges[that_net][:end_name] + '_clustered.tsv'
	#network_edges.to_csv(folder + '/' + cled_name, sep = '\t')

# 	plot(communities[1])
# 	show()
	
	
	
	
#	print(network_mat.values[:,1:].shape)
	
#	net_graph = nx.from_numpy_matrix(network_mat.values[:,1:])

# communities = com_clust(network_mat)
# print(communities)

#print(var(network_mat.values[:,1:]))


##### Plots of stats (rough) ###################
# xpos = arange(len(sample_types))
# fig, (ax1,ax2) = plt.subplots(nrows = 2)
# plt.xticks(rotation=90)
# 
# ax1.bar(xpos, mean_weights, align='center')
# ax2.bar(xpos, edge_props, align='center')
# 
# ax1.set_xticks(xpos)
# ax1.set_xticklabels(sample_types)
# ax1.set_title('Mean Edge Weights')
# 
# ax2.set_xticks(xpos)
# ax2.set_xticklabels(sample_types)
# ax2.set_title('Propotion of Edges')
# 
# plt.show()



##############################################################################
#######										##################################
######		To Do 							##################################
#######										##################################
##############################################################################
#
#
# - Implement network clustering and compare to grouping by sample type.













