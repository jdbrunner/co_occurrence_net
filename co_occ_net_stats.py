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

###### Conductance of cuts ##################
##
# - Calculate conductance of cuts that cut out a type or set of types.
#	- This will give an idea about the strength of intertype nodes - because they will
#		always be the ones cut. Also allows us to group types (i.e. gums with cheek)

##### Plots of stats (rough) ###################
xpos = arange(len(sample_types))
fig, (ax1,ax2) = plt.subplots(nrows = 2)
plt.xticks(rotation=90)

ax1.bar(xpos, mean_weights, align='center')
ax2.bar(xpos, edge_props, align='center')

ax1.set_xticks(xpos)
ax1.set_xticklabels(sample_types)
ax1.set_title('Mean Edge Weights')

ax2.set_xticks(xpos)
ax2.set_xticklabels(sample_types)
ax2.set_title('Propotion of Edges')

plt.show()



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













