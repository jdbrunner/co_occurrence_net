#########################################################################################
#																						#		
#				Co-Occurrence Network Builder											#
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
import itertools as iter
from scipy.special import binom as choose
from scipy.stats import binom
import re


#####Variables
##User input, should be entered as an argument
csv_name = sys.argv[1]
level = sys.argv[2]
net_name = sys.argv[3]
sample_types = sys.argv[4]
if sample_types == 'True':
	sample_types = True
else:
	sample_types = False

############# Functions
#Function that counts the number of times both r1 and r2 are within a range, and then 
#returns the fraction of times this occurs. The range is half open: (a,b]
def both_occ(r1,r2, lbd = 0, ubd = 1):
	co_occ = count_nonzero(((r1>lbd) & (r1<=ubd) & (r2>lbd) & (r2<=ubd)))
	co_occ = float(co_occ)
	return co_occ

#Function that counts the number of times the two occur in the same abundance range.
#rel = true does this relative to the maximum abundance
def both_same_bin(r1,r2, lthresh, numthresh, rel = True):
	threshes = linspace(0, 1, numthresh)
	threshes[0] = lthresh
	k = len(threshes)
	if rel:
		r1m = max(r1)
		r2m = max(r2)
		r1 = r1/r1m
		r2 = r2/r2m	
	both_in_bin = [both_occ(r1,r2,lbd = threshes[i],ubd = threshes[i+1]) for i in range(k-1)]
	tot_same_bin = sum(both_in_bin)
#	same_bin_frac = float(tot_same_bin)/float(len(r1))
	return tot_same_bin
	

#Function that classifies nodes by which type of sample they have the highest abundance in.
#I need to consider how to deal with nodes that appear at high levels in more than one 
#type of sample.
def color_picker(r):
	numcols = len(r)
	most_found = argmax(r)
	strgth = r[most_found] - mean(r)
	std_dev = std(r)
	if strgth >= 1.5*std_dev:
		colors = linspace(0,255,numcols)
		numwhere = where(r.index == argmax(r))
		color_picked = colors[numwhere]
		rgb_value = cm.rainbow(color_picked.astype(int))
		hex_value = matplotlib.colors.rgb2hex(rgb_value[0])
		return [most_found, hex_value]
	else:
		grey = matplotlib.colors.rgb2hex([0.80000001192092896, 0.80000001192092896, 0.80000001192092896, 1.0])
		return ['Mixed', grey]

def matchyn(a,b):
	if a==b:
		return a
	else:
		grey = matplotlib.colors.rgb2hex([0.80000001192092896, 0.80000001192092896, 0.80000001192092896, 1.0])
		return ['Intertype', grey]
		
def occ_probs(abund_array,lthresh, numthresh, rel = True):
	'''Calculate probability of occurrence at an abundance level in a random graph, binomial distribution
	with parameters edge degree and sample degree/total edges, as in \cite{coocc}.'''
	#Make the bipartite graph by creating a list for each sample (converting 
	#to discrete values based on bin number)
	if rel:
		for rw in abund_array.index:
			maxab = max(abund_array.loc[rw][1:])
			abund_array.loc[rw,1:] = abund_array.loc[rw][1:]/maxab
	samp = abund_array.columns[1:]
	bins = ['_b'+ str(i) for i in range(1,numthresh+1)]
	bigraph = pd.DataFrame(abund_array['TAXA'],columns = ['TAXA'])
	threshes = linspace(0, 1, numthresh)
	threshes[0] = lthresh
	for this_samp in samp:
		for j in range(1,numthresh):
			bigraph[this_samp+bins[j-1]] = ((abund_array[this_samp] > threshes[j-1]) & (abund_array[this_samp] <= threshes[j])).astype(int)
	samp_degs = sum(bigraph,0)[1:]
	org_degs = sum(bigraph,1)
	total_edges = sum(samp_degs)
	bin_prob = samp_degs/total_edges
	occ_prob = array([[1 - (1-p)**n for p in bin_prob] for n in org_degs])
	return occ_prob

############## NOT FOR USE WITHOUT A COMPUTER FROM 2117	
def random_coocc_prob(occ,wij,i,j):
	'''Calculate a poisson-binomial...P(X > wij) where X is the number of times
	i and j co occur in random graph (X is a random variable)'''
	#this is crazy slow. Luckily theres a paper with an idea about approximating it!
	#occ is a matrix with entry il the probability i occurs in sample l under the 
	#random graph
	coocc = occ[i]*occ[j] #vector with p_il*p_jl
	no_coocc = 1 - coocc
	numsamp = len(occ[0])
	terms = range(wij,numsamp)
	prob = 0
	for k in terms:
		choices = array(list(iter.combinations(range(numsamp),k)))
		for ch in choices:
			unch = delete(range(numsamp),ch)
			prob = prob + prod(coocc[ch])*prod(no_coocc[unch])
	return prob
##########################################

def approx_rand_prob(occ,wij,i,j):
	'''Calculate an approximation of a poisson-binomial...P(X > wij) where X 
	is the number of timesi and j co occur in random graph (X is a random variable)
	The paper \cite{coocc} calls it bi-binomial because its as if you had two 
	probabilities in the above'''
	coocc = occ[i]*occ[j] #vector with p_il*p_jl
	N = len(coocc)
	mu = sum(coocc)
	sig2 = mu - sum(coocc**2)
	pa = mu/N
	N2 = round(mu**2/(mu - sig2))
	N1 = N-N2
	S2 = (mu-sig2)/N - (mu/N)**2
	p1 = pa - sqrt((N2*S2/N1))
	p2 = pa + sqrt((N1*S2/N2))
	terms = range(int(wij),10)
	prob = 0
	###for k = wi1,...,N, calculate P(X = k) and add it on
	for k in terms:
		ztok = array(range(k+1))
		kminj = k - ztok
		n1side = binom.pmf(ztok, N1, p1)
		n2side = binom.pmf(kminj, N2, p2)
		prob = prob + sum(n1side*n2side)
	return prob

######Import the abundance matrix
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 
abundance_array_full = pd.read_csv(csv_name, sep = ' ')
#print(abundance_array)
#Let's remove taxa that aren't present at all.
for ind in abundance_array_full.index:
	if max(abundance_array_full.loc[ind][2:]) == 0:
		abundance_array_full = abundance_array_full.drop(ind)

abundance_array_full.index = range(len(abundance_array_full))		

#We want to make a network at a specific level, or set of levels
if level == 'all':
	level = array(abundance_array_full['LEVEL'].values)
	level = array(unique(level))
else:
	level = array([level])

#need to check that the levels given as argument actually exist
levels_allowed = array(abundance_array_full['LEVEL'].values)
levels_allowed = array(unique(levels_allowed))

for i in level:
	print(i)
	if not i in levels_allowed:
		print(levels_allowed)
		sys.exit('Bad level argument. Choose a level from the list above or type all')

#Combine samples of the same kind (that is, the same part of the body) so that we can 
#color according to where abundance of a genome is highest.
diff_samps_types = unique([name[:-10] for name in array(abundance_array_full.keys())[2:]])
data_samps_types = transpose([sum([abundance_array_full[samp] for samp in abundance_array_full.keys()[2:]
							 if smp_type in samp],axis = 0) for smp_type in diff_samps_types])

samp_type_abund = pd.DataFrame(data_samps_types, columns = diff_samps_types)
samp_type_abund['TAXA'] = abundance_array_full['TAXA']

the_levels = abundance_array_full['LEVEL'].values
the_taxa = abundance_array_full['TAXA'].values

ab_arrays = []	
if sample_types:
	for samp in diff_samps_types:
		slist = []
		ab_samp = pd.DataFrame( the_levels, columns = ['LEVEL'], index = abundance_array_full.index)
		ab_samp['TAXA'] = the_taxa
		selecter = abundance_array_full.columns.map(lambda x: bool(re.search(samp,x)))
		ab_samp[abundance_array_full.columns[where(selecter)]] = abundance_array_full[abundance_array_full.columns[where(selecter)]]
		ab_arrays = ab_arrays + [ab_samp]
else:
	ab_arrays = [abundance_array_full]

adjacency_frames = dict()
source_target_frames = dict()
adjecency_frames_pre = dict()
source_target_frames_pre = dict()
#colors = dict()
for i in level:
	for ii in range(len(ab_arrays)):
		if sample_types:
			stype = diff_samps_types[ii]
		else:
			stype = 'all'
		abundance_array = ab_arrays[ii]
		#start by getting the indices of the members of the level
		in_level = where(abundance_array['LEVEL'] == i)[0]
		#then the part of the dataframe containing just that level
		lvl_abundance_array = abundance_array.iloc[in_level]
		lvl_abundance_array = lvl_abundance_array.drop('LEVEL',axis = 1)
		#create np array of abundance values (dropping labels) within the level
		ab_np_array1 = lvl_abundance_array.values[:,1:]
		not_seen = where(sum(ab_np_array1,1) == 0)[0]
		n_seen_ind = lvl_abundance_array.index[not_seen]
		in_level = delete(in_level,not_seen)
		#prune off unseen organisms
		lvl_abundance_array = lvl_abundance_array.drop(n_seen_ind)
		ab_np_array = lvl_abundance_array.values[:,1:]	
		occur_probs = occ_probs(lvl_abundance_array,0.05,5)
	
	
		#Create numpy array of co-occurrence fractions. This is the incidence matrix for our weighted
		#graph.
		adj_matrix_pre = asarray([[0 if x == y else both_same_bin(ab_np_array[x],ab_np_array[y], 0.05, 5) 
									for x in range(len(ab_np_array))] for y in range(len(ab_np_array))])
		degs = sum(adj_matrix_pre,1)
		unconnected = where(degs == 0)
		adj_matrix_pre = delete(adj_matrix_pre,unconnected,0)
		adj_matrix_pre = delete(adj_matrix_pre,unconnected,1)
		in_level = delete(in_level,unconnected)
		
		adjecency_frames_pre[i + '_' + stype] = pd.DataFrame(adj_matrix_pre, index = abundance_array['TAXA'][in_level], columns = abundance_array['TAXA'][in_level])

	
		stringency = 0.05
	
		N = len(ab_np_array[0])
		tot_seen = [sum([abd != 0 for abd in row]) for row in ab_np_array]
		tot_seen = delete(tot_seen,unconnected)
	
		adj_size = adj_matrix_pre.shape
		adj_matrix = zeros(adj_size)
		for k in range(adj_size[0]):
			for j in range(adj_size[1]):
				if adj_matrix_pre[k,j] != 0:
					p = approx_rand_prob(occur_probs,adj_matrix_pre[k,j],k,j)
					if p <= stringency:
						adj_matrix[k,j] = adj_matrix_pre[k,j] 
	
		degs2 = sum(adj_matrix,1)
		unconnected2 = where(degs2 == 0)
		adj_matrix = delete(adj_matrix,unconnected2,0)
		adj_matrix = delete(adj_matrix,unconnected2,1)
		tot_seen_2 = delete(tot_seen,unconnected2)
				
		in_level_2 = delete(in_level,unconnected2)
		
		adjacency_frames[i + '_' + stype] = pd.DataFrame(adj_matrix, index = abundance_array['TAXA'][in_level_2], columns = abundance_array['TAXA'][in_level_2])
	
		toedge = True
		if toedge:
			#create a list of edges - source in one column, target in the other. This is an alternative to the adjacency matrix
			#that might make it easier to load in to cytoscape.
			num_edge = count_nonzero(adj_matrix)
	# 		print(i)
	# 		print(num_edge/2)
			source_target_data = []#empty([num_edge,3], dtype = 'string')	
			for l in range(len(adj_matrix)):
				for k in range(l+1):
						if adj_matrix[l,k] != 0:
							s_type1 = color_picker(samp_type_abund.iloc[in_level_2[l]][:-1])
							s_type2 = color_picker(samp_type_abund.iloc[in_level_2[k]][:-1])
							edge_samp = matchyn(s_type1,s_type2)
							#include the "reverse" edge as well so that cytoscape doesn't have gaps in 
							#it's classification of nodes.
							edge =  [abundance_array['TAXA'][in_level_2[l]], abundance_array['TAXA'][in_level_2[k]],
														str(adj_matrix[l,k]), str(tot_seen[l]),str(tot_seen[k]),
														s_type1[0],s_type2[0],edge_samp[0],s_type1[1],s_type2[1],edge_samp[1],N]
							rev_edge = [abundance_array['TAXA'][in_level_2[k]], abundance_array['TAXA'][in_level_2[l]],
														str(adj_matrix[l,k]),str(tot_seen_2[k]),str(tot_seen_2[l]),
														s_type2[0],s_type1[0],edge_samp[0],s_type2[1],s_type1[1],edge_samp[1],N]
							source_target_data += [edge]
							source_target_data += [rev_edge]
		
			source_target_data_pre = []#empty([num_edge,3], dtype = 'string')	
			for l in range(len(adj_matrix_pre)):
				for k in range(l+1):
						if adj_matrix_pre[l,k] != 0:
							s_type1 = color_picker(samp_type_abund.iloc[in_level[l]][:-1])
							s_type2 = color_picker(samp_type_abund.iloc[in_level[k]][:-1])
							edge_samp = matchyn(s_type1,s_type2)
							#include the "reverse" edge as well so that cytoscape doesn't have gaps in 
							#it's classification of nodes.
							edge =  [abundance_array['TAXA'][in_level[l]], abundance_array['TAXA'][in_level[k]],
														str(adj_matrix_pre[l,k]),str(tot_seen[l]),str(tot_seen[k]),s_type1[0],s_type2[0],edge_samp[0],s_type1[1],s_type2[1],edge_samp[1]]
							rev_edge = [abundance_array['TAXA'][in_level[k]], abundance_array['TAXA'][in_level[l]],
														str(adj_matrix_pre[l,k]),str(tot_seen[k]),str(tot_seen[l]),s_type2[0],s_type1[0],edge_samp[0],s_type2[1],s_type1[1],edge_samp[1]]
							source_target_data_pre += [edge]
							source_target_data_pre += [rev_edge]
				
			#turn list into a dataframe
			source_target_frames[i + '_' + stype] = pd.DataFrame(source_target_data, columns = ['source','target','weight','source_freq',
								'target_freq','first_sample','second_sample','edge_sample','fscolor','sscolor','edcolor','num_samps'])
			source_target_frames_pre[i + '_' + stype] = pd.DataFrame(source_target_data_pre, columns = ['source','target','weight','source_freq',
								'target_freq','first_sample','second_sample','edge_sample','fscolor','sscolor','edcolor'])
	
		#Assign a color to each node based on where they are most abundant (what sample type)

#Save adjacency_frames to save the (weighted) adjacency matrix. This can be loaded into 
#cytoscape. Save source_target_frames to save a list of edges. This can also be loaded into
#cytoscape and provides a way to include edge or node attributes.
save = True
if save:
	for j in source_target_frames_pre.keys():
		flname1 = net_name+'_'+j+'_coin_adj.tsv'
		flname2 = net_name+'_'+j+'_coin_list.tsv'
		flname3 = net_name+'_'+j+'_coocc_adj.tsv'
		flname4 = net_name+'_'+j+'_coocc_list.tsv'
		adjecency_frames_pre[j].to_csv(flname1, sep = '\t')
		source_target_frames_pre[j].to_csv(flname2, sep = '\t')
		adjacency_frames[j].to_csv(flname3, sep = '\t')
		source_target_frames[j].to_csv(flname4, sep = '\t')


#### TO DO #####################
#
# - Decide how to classify nodes that appear in multiple types of sample at 
#				high frequency
#
# - Change co-occurrence network to account for randomness
#















