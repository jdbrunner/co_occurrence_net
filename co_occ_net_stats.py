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
	
def exp_cut_edges(network,types, p = 0.5):
	#first grab all the edges within the relavent subnetwork
	subedges = network.iloc[where([type in types for type in network['edge_sample']])]
	#Now count the nodes in the subnetwork
	num_nodes = len(unique(subedges['source']))
	num_edges = len(subedges)/2
	#Now count the nodes not in the subnetwork
	tot_nodes = len(unique(network['source'])) - num_nodes
	#number of possible intertype edges
	all_inter = num_nodes*tot_nodes
	expected = p*all_inter
	actual = len(where([type == 'Intertype' for type in network['edge_sample']]))
	return [expected, actual]
	
these_types = sys.argv[2]
print exp_cut_edges(network_edges,these_types)

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













