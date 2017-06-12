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
from scipy import misc

#####Variables
##User input, should be entered as an argument
csv_name = sys.argv[1]


######Import the network as list of edges
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 
network_edges = pd.read_csv(csv_name, sep = '\t')
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
	
###### Number of Edges vs Random ############
##
# - Compare number of types of edges vs expected number of that type in a random graph
#   - Use ER random graph
#	- For edges within a sample type this just number of edges vs E(number of edges in 
# 		an ER random graph with number of nodes = number of nodes of that sample type)
#	- For intertype edges this will be harder. I'll have to check this but it should be
#		E(edges of ER random graph) - \sum_{types} E(edges of subgraph)
# 
these_types = sys.argv[2]
if ',' in these_types:
	commas = where([letter == ',' for letter in these_types])
	commas = append(commas,len(these_types)-1)
	commas = insert(commas,0,0)
	print commas
	the_types = []
	for i in range(1,len(commas)):
		the_types = the_types + ["".join(these_types[commas[i-1]+1:commas[i]])]
	these_types = the_types
else:
	these_types = ["".join(these_types[1:-1])]


def edge_prob(network):
	#calculate probability of seeing an edge in a graph with that many edges
	#if edge prob is p and we see m edges, p \approx  2m/n^2 (n nodes)
	n = len(unique(network['source']))
	m = len(network)/2
	p = float(2*m)/float(n**2 - n)
	return p
	
e_prob = edge_prob(network_edges)


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
	
print random_sub_graph(network_edges,these_types, p = e_prob)
	
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
	
print exp_cut_edges(network_edges,these_types, p = e_prob,  between = False)

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
# - Compare edges to random graph - decide on type of random graph (ie edge distribution)
#
# - Calculate conductance from sample to type to rest of graph:
#		By calculating conductance of cutting out that sample type, weighted by number
#			of nodes of that sample.
#
# - Implement network clustering and compare to grouping by sample type.













