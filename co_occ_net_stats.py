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
#import networkx as nx

#####Variables
##User input, should be entered as an argument
folder = sys.argv[1]
flder = os.listdir(folder)
coin_matrices = [flder[j] for j in where([i[-12:] == 'coin_adj.tsv' for i in flder])[0]]
coin_edges = [flder[j] for j in where([i[-13:] == 'coin_list.tsv' for i in flder])[0]]
coocc_matrices = [flder[j] for j in where([i[-13:] == 'coocc_adj.tsv' for i in flder])[0]]
coocc_edges = [flder[j] for j in where([i[-14:] == 'coocc_list.tsv' for i in flder])[0]]
numnets = len(coin_matrices)

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


def edge_prob(network):
	#calculate probability of seeing an edge in a graph with that many edges
	#if edge prob is p and we see m edges, p \approx  2m/n^2 (n nodes)
	n = len(unique(network['source']))
	m = len(network)/2
	p = float(2*m)/float(n**2 - n)
	return p
	

#count the nodes of a type
def nodes_in_sub(network,type):
	nodes = network.iloc[where(network['first_sample']==type)]
	num_nodes = len(unique(nodes['source']))
	return num_nodes

### Calculate probability that the subnetwork has as many edges as it has
def random_sub_graph(network,types, p = 0.5):
	#first grab all the edges within the relavent subnetwork
	subedges = network.iloc[where([type in types for type in network['edge_sample']])]
	#Now count the nodes in the subnetwork
	num_nodes = len(unique(subedges['source']))
	num_edges = len(subedges)/2
	#at most num_nodes choose 2 edges
	full_graph = (num_nodes*(num_nodes - 1)/(float(2)))
	expected = p*full_graph
	#prob = float(misc.comb(full_graph,num_edges))*(p**num_edges)*((1-p)**(full_graph - num_edges))
	#diff = expected - num_edges
	return [expected, num_edges]
		
### Calculate the probability of seeing that many intertype edges, or edges between some
#		certain set of types, or out of some set of types.
def exp_cut_edges(network, types, p = 0.5, between = True):
	if between == False and types == 'all':
	#you can't count the edges out of all the types
		return [0,0]
	else:
		#find the edges between types
		inter_edges_loc = where(network['edge_sample'] == 'Intertype')
		edges = network.iloc[inter_edges_loc]
		if types == 'all':
		#need the list of types
			types = delete(unique(network['first_sample'].values),where(network['first_sample'].values == 'Intertype'))
		else:
			#narrow down if you want to
			if between:
				betwn_loc_1 = where([type in types for type in edges['first_sample']])
				betwn_loc_2 = where([type in types for type in edges['second_sample'].iloc[betwn_loc_1]])
				betwn_loc = betwn_loc_1[0][betwn_loc_2[0]]
				edges = edges.iloc[betwn_loc]
			else:
				out_loc_1 = where([type in types for type in edges['first_sample']])
				out_loc_2 = where([type not in types for type in edges['second_sample'].iloc[out_loc_1]])
				out_loc = out_loc_1[0][out_loc_2[0]]
				edges = edges.iloc[out_loc]
		if between:
			#if we are taking edges between, we have double counted
			num_edges = len(edges)/2
			#number of nodes of each type
			nodes_per = [nodes_in_sub(network,type) for type in types]
			#then compute the total pairs of nodes between
			int_poss = 0
			for i in range(len(nodes_per)):
				for j in range(i,len(nodes_per)):
					int_poss = int_poss + nodes_per[i]*nodes_per[j]
			#and multiply by edge probability
			int_expected = int_poss*p
		else:
			num_edges = len(edges)
			#number of nodes of each type
			in_nodes = 0
			for type in types:
				in_nodes = in_nodes + nodes_in_sub(network,type)
			#and number of nodes of none of these types
			#all the types
			out_nodes = len(unique(network['source'])) - in_nodes
			out_tot = in_nodes*out_nodes
			int_expected = out_tot*p
		return [int_expected, num_edges]
	

###### Conductance of cuts ##################
##
# - Calculate conductance of cuts that cut out a type or set of types.
#	- This will give an idea about the strength of intertype nodes - because they will
#		always be the ones cut. Also allows us to group types (i.e. gums with cheek)

def cut_cond(network,types):
	###Isolate relevant edges (get their integer index)
	sub_edges = where([type in types for type in network['first_sample']])
	##Pick out edges to be cut - this is the index in rel_edges of their index in network (woof) 
	cut_edges = where([type not in types for type in network['second_sample'].iloc[sub_edges]])
	edges_to_be_cut = network.iloc[sub_edges[0][cut_edges]]
	#
	sub_only_edges_loc =  where([type in types for type in network['second_sample'].iloc[sub_edges]])
	sub_only_edges = network.iloc[sub_edges[0][sub_only_edges_loc]]
	#
	non_sub_edges_pre =  where([type not in types for type in network['first_sample']])
	non_sub_edges_loc =  where([type not in types for type in network['second_sample'].iloc[non_sub_edges_pre]])
	non_sub_edges = network.iloc[non_sub_edges_pre[0][non_sub_edges_loc]]
	#
	cut_cond_num = sum(edges_to_be_cut['weight'])
	a_non_sub = sum(non_sub_edges['weight'])
	a_sub = sum(sub_only_edges['weight'])
	cut_cond_den = min(a_non_sub, a_sub)
	if cut_cond_den != 0:
		return float(cut_cond_num)/float(cut_cond_den)
	else:
		return 'Empty Cut'
	
	
############# Clustering #######################################
#####
### - Cluster nodes using community detection cluster and 
#		spectral clustering (I would like to compare)





def com_clust(network):
	#Clustering using Girvan and Newman algortithm. This means we try to maximize the 
	#quantity \emph{modularity} over all possible community groupings. Modularity is
	#the quantity: (Fraction of edges that are inside a community) - (expected fraction in a random graph)
	#I'd like to weight this, so I'll maximize
	#((\delta is Kronecher, d is sum of weights of edges on that vertex (generalized degree))
	#Also, the undirectedness means I only have to take the sum over the subdiagonal. 
	adjmat = network.values[:,1:]
	sum_weight = sum(adjmat)
	ai = 1/(sum_weight)*array([sum(row) for row in adjmat])
#	Qterms = 1/(sum_weight)*adjmat - outer(ai,ai)
	deltaQ = 1/(sum_weight)*adjmat - outer(ai,ai)#Qterms
#	oneclus = sum(deltaQ)
	maxDQ = unravel_index(deltaQ.argmax(), deltaQ.shape)
#	removed = []
	Q = [sum(diag(deltaQ))]
	cluster_list = []
	for i in range(adj_mat.shape[0]):
		cluster_list = cluster_list+[[i]]
	while deltaQ[maxDQ] > 0:#deltaQ.shape[0] > 1:
#		removed = removed+[maxDQ]
		cluster_list[maxDQ[0]] = cluster_list[maxDQ[0]] + cluster_list[maxDQ[1]]
		cluster_list.remove(cluster_list[maxDQ[1]])
		Q = Q + [Q[-1] + 2*deltaQ[maxDQ]]
		deltaQ[maxDQ[0]] = deltaQ[maxDQ[0]]+deltaQ[maxDQ[1]]
		deltaQ[:,maxDQ[0]] = deltaQ[:,maxDQ[0]]+deltaQ[:,maxDQ[1]]
		deltaQ[maxDQ[0],maxDQ[0]] = deltaQ[maxDQ[0],maxDQ[0]]
		deltaQ = delete(deltaQ, maxDQ[1],0)
		deltaQ = delete(deltaQ, maxDQ[1],1)
		mask = array(deltaQ)
		fill_diagonal(mask, -inf)
		maxDQ = unravel_index(mask.argmax(), deltaQ.shape)
# 		if deltaQ[maxDQ] <= 0:
# 			print(deltaQ.shape)
	#cluster_list list of lists - each sublist is a cluster. the sublists contain the index
	#of the nodes in that cluster. It will be nicer to have a list where entry i tells us 
	#what cluster node i is in.
	comm_clusts = [where([node in cluster for cluster in cluster_list])[0][0] for node in range(adjmat.shape[0])]
	return [comm_clusts,Q]
# 		
		
# 	
def spectral_cluster(adj_mat):
	'''Clustering using spectral clustering.'''
	#construct the graph laplacian. Need degree matrix first
	degs = sum(adj_mat,1)
	D = diag(degs)
	Din = diag(1/degs)
	L = D - adj_mat
	Lrw = dot(Din,L)
	[evals,evects] = eigh(Lrw.astype(float))
	most_gaps = min([len(evals),20])
	sgaps = [evals[i] - evals[i-1] for i in range(1,most_gaps)]
	num_clus = 5 + argmax(sgaps[5:])
	V = evects[:,:num_clus]
	Vwhite = cluster.vq.whiten(V)
	cents = cluster.vq.kmeans(Vwhite,num_clus)[0]
	spect_clusts = cluster.vq.vq(Vwhite, cents)
	#returns list - entry i tells us what cluster node i is in
	return spect_clusts[0]




def color_picker2(r):
	'''Takes in vector r and maps to hex colors'''
	muncols = len(r)
	colors = linspace(0,255,muncols)
	collist = []
	for i in range(muncols):
		colori = colors[i]
		rgb_value = cm.rainbow(colori.astype(int))
		hex_value = matplotlib.colors.rgb2hex(rgb_value)
		collist = collist + [hex_value]
	return collist



		
######Import the network as list of edges
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 

for this_net in range(numnets):
	network_edges = pd.read_csv(folder + '/' + coin_edges[this_net], sep = '\t')
	network_edges.index = network_edges['source']
	network_mat = pd.read_csv(folder + '/' + coin_matrices[this_net],sep = "\t")

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
	cled_name = coin_edges[this_net][:end_name] + '_clustered.tsv'
	network_edges.to_csv(folder + '/' + cled_name, sep = '\t')


for that_net in range(numnets):
	network_edges = pd.read_csv(folder + '/' + coocc_edges[that_net], sep = '\t')
	network_edges.index = network_edges['source']
	network_mat = pd.read_csv(folder + '/' + coocc_matrices[that_net],sep = "\t")

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
	cled_name = coocc_edges[that_net][:end_name] + '_clustered.tsv'
	network_edges.to_csv(folder + '/' + cled_name, sep = '\t')

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













