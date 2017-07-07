#########################################################################################
#																						#		
#				Clean up a network file for import into cytoscape						#
#																						#
#										Jim Brunner 									#
#																						#
#																						#
#																						#
#########################################################################################

import pandas as pd
from pylab import *
import sys

net_name = sys.argv[1]
net = pd.read_csv(net_name, sep = '\t')

for col in unique(net.columns.values):
	if 'Unnamed' in col:
		net.drop(col,axis = 1, inplace = True)



to_drp1 = ['source.1','num_samps']
for coll in to_drp1:
	if coll in net.columns:
		net.drop(coll,axis = 1, inplace = True)

net.to_csv(net_name,sep = '\t')

to_drp2 = ['target_freq','sscolor','second_sample']
for colll in to_drp2:
	if colll in net.columns:
		net.drop(colll,axis = 1, inplace = True)


edge_info = net[['source','target','weight','edge_sample','edcolor']]

node_info = net.drop(['target','weight','edge_sample','edcolor'], axis = 1)


slsh = where([let == '/' for let in net_name])[0][-1]

edge_info.to_csv(net_name[:slsh]+'/cyto_input/'+net_name[slsh+1:-4]+'_edges.tsv', sep = '\t')
node_info.to_csv(net_name[:slsh]+'/cyto_input/'+net_name[slsh+1:-4]+'_node_atts.tsv', sep = '\t')



