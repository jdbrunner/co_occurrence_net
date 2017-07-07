# co_occurrence_net
Co-occurrence network builder and analysis based on this network

co_occurrence.py is a python script to produce a co-occurrence network that can be given to Cytoscape. Outputs a .tsv file with list of edges. Input is data in merged_assignment.txt file.

co_occ_net_stats.py does some network statistics and clustering

color_key.py just plots colors

rand_samp.py produces a "sample", both random and taken from a column of merged_assignment.txt

sample_analysis.py takes a sample and a network and adds columns to the network file corresponding to that sample - inclusing "likelihood of presence rankings" from diffusion based analysis

cleanup.py takes a network file (list of edges and node attributes) and returns a file with just the source, target, and edge attributes, as well as a second file with source attributes.

To import a network into cytoscape, run cleanup.py on it, then import the edges file (import network from file). Next, import the nodes file (import table from file) as a node attribute table.

There are also two bash scripts which will do these things. The first makes and clusters networks. The second creates and analyzes samples.