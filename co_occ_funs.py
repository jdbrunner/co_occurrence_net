#########################################################################################
#																						#		
#				Module for co-occurrence network analysis								#
#																						#
#										Jim Brunner 									#
#																						#
#																						#
#																						#
#########################################################################################





from pylab import *
import pandas as pd
import itertools as iter
from scipy.special import binom as choose
from scipy.stats import binom
from numpy.random import binomial as bino
import re
import numpy.ma as mask
from scipy import misc, sparse, cluster
from random import sample




####################################  Building a co-occurrence network			 #########
###########################################################################################
######################

def both_occ(r1,r2, lbd = 0, ubd = 1):
	'''Function that counts the number of times both r1 and r2 are within a range, and then 
	returns the fraction of times this occurs. The range is half open: (a,b]'''
	co_occ = count_nonzero(((r1>lbd) & (r1<=ubd) & (r2>lbd) & (r2<=ubd)))
	co_occ = float(co_occ)
	return co_occ


def both_same_bin(r1,r2, lthresh, numthresh, rel = True):
	'''Function that counts the number of times the two occur in the same abundance range.
	rel = true does this relative to the maximum abundance'''
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


#I need to consider how to deal with nodes that appear at high levels in more than one 
#type of sample.
def color_picker(r):
	'''Function that classifies nodes by which type of sample they have the highest abundance in.'''
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


def mc_dot_prod(N,P,W,num_samps = 1000):
	'''MC approximation for P(X\cdot Y > w) where X is a random vector of binomial(n1,p1)
	and Y is a random vector of binomial(n2,p2), where p1,p2 are vectors'''
	mc_samples = zeros((num_samps,N.shape[0],N.shape[0]))
	N = N.astype(int)
	for s in range(num_samps):
		the_rand_mat = bino(N,P)
		rsumd = sum(the_rand_mat,1)
		rsumd[where(rsumd == 0)] = 1
		normed_rand_m = dot(diag(1/rsumd),the_rand_mat)
		rand_weights = dot(normed_rand_m,transpose(normed_rand_m))
		comp_mat = rand_weights - W
		mc_samples[s][where(comp_mat >= 0)] =  1
	return mean(mc_samples,axis = 0)#,var(mc_samples,axis = 0).max()]
	
def mc_pearson(N,P,W,samp_types = None,data_header = None,num_samps = 1000):
	'''MC approximation for pearson coefficient where X is a random vector of binomial(n1,p1)
	and Y is a random vector of binomial(n2,p2), where p1,p2 are vectors. Expected value is identity.'''
	mc_samples = zeros((num_samps,N.shape[0],N.shape[0]))
	N = N.astype(int)
	for s in range(num_samps):
		the_rand_mat = bino(N,P)
		if samp_types != None:
			for idx in samp_types:
				selec = data_header.map(lambda x: bool(re.search(idx,x)))
				selllec = where(selec[1:])
				meanss = mean(the_rand_mat[:,selllec[0]],axis = 1)
				vars = std(the_rand_mat[:,selllec[0]],axis = 1)
				the_rand_mat[:,selllec[0]] = the_rand_mat[:,selllec[0]] - outer(meanss,ones(len(selllec[0])))
				vars[where(vars == 0)] = 1
				for idx2 in range(len(the_rand_mat)):
					the_rand_mat[idx2,selllec[0]] = the_rand_mat[idx2,selllec[0]]/vars[idx2]
				rand_weights = (1/the_rand_mat.shape[1])*dot(the_rand_mat,transpose(the_rand_mat))
		else:
			meanss = mean(the_rand_mat,axis = 1)
			vars = std(the_rand_mat,axis =1)
			vars[where(vars == 0)] = 1
			the_rand_mat = dot(diag(vars),transpose(transpose(the_rand_mat) - meanss))
			rand_weights = (1/the_rand_mat.shape[1])*dot(the_rand_mat,transpose(the_rand_mat))
# 		print(rand_weights.shape)
# 		print(W.shape)
		comp_mat = rand_weights - W #- eye(rand_weights.shape[0])
		mc_samples[s][where(comp_mat >= 0)] =  1
	return mean(mc_samples,axis = 0)#,var(mc_samples,axis = 0).max()]
	
def min_nz(arr, rows = False):
	'''Find minimum non-zero value'''
	mask_arr = mask.masked_values(arr, 0, copy = False)
	if rows:
		return mask_arr.min(axis = 1).data
	else:
		return mask_arr.min()
		
		
#########################################
####################################  Analysis of a whole network			 #########
###########################################################################################
######################

def edge_prob(network):
	'''calculate probability of seeing an edge in a graph with that many edges
	if edge prob is p and we see m edges, p \approx  2m/n^2 (n nodes)'''
	n = len(unique(network['source']))
	m = len(network)/2
	p = float(2*m)/float(n**2 - n)
	return p
	

#
def nodes_in_sub(network,type):
	'''count the nodes of a type'''
	nodes = network.iloc[where(network['first_sample']==type)]
	num_nodes = len(unique(nodes['source']))
	return num_nodes

###
def random_sub_graph(network,types, p = 0.5):
	''' Calculate probability that the subnetwork has as many edges as it has'''
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
		
### 
def exp_cut_edges(network, types, p = 0.5, between = True):
	'''Calculate the probability of seeing that many intertype edges, or edges between some
	certain set of types, or out of some set of types.'''
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
	'''Calculate conductance of cuts that cut out a type or set of types'''
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
	'''Clustering using Girvan and Newman algortithm. This means we try to maximize the 
	quantity \emph{modularity} over all possible community groupings. Modularity is
	the quantity: (Fraction of edges that are inside a community) - (expected fraction in a random graph)
	I'd like to weight this, so I'll maximize
	((\delta is Kronecher, d is sum of weights of edges on that vertex (generalized degree))
	Also, the undirectedness means I only have to take the sum over the subdiagonal. '''
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
	for i in range(adjmat.shape[0]):
		cluster_list = cluster_list+[[i]]
	while deltaQ[maxDQ] > 10**(-10):#deltaQ.shape[0] > 1:
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
	[evals,evects] = eig(Lrw.astype(float))
	ord = argsort(evals)
	evals = evals[ord].real
	evects = evects[:,ord].real
	max_cl = len(argwhere(evals.round(10) == 0)) + 10
	most_gaps = min([len(evals),max_cl])
	sgaps = [evals.round(10)[i] - evals.round(10)[i-1] for i in range(1,most_gaps)]
	num_clus = 5 + argmax(sgaps[5:])
	V = evects[:,:num_clus].astype(float)
	Vwhite = cluster.vq.whiten(V)
	cents = cluster.vq.kmeans(Vwhite,num_clus)[0]
	spect_clusts = cluster.vq.vq(Vwhite, cents)
	#returns list - entry i tells us what cluster node i is in
	return spect_clusts[0]




def color_picker2(r,the_map = cm.rainbow, weighted = False):
	'''Takes in vector r and maps to hex colors'''
	muncols = len(r)
	if weighted:
		ma = max(r)
		mi = min(r)
		if ma - mi != 0:
			slp = 255/(ma - mi)
			colors = (r*slp - slp*mi).astype(int)
		else:
			colors = zeros(muncols,dtype = int)
	else:
		colors = linspace(0,255,muncols).astype(int)
	rgb_value_array = the_map(colors)
	collist = []
	for i in range(muncols):
		rgb_value = rgb_value_array[i]
		hex_value = matplotlib.colors.rgb2hex(rgb_value)
		collist = collist + [hex_value]
	return collist




####################################  Determining likelihood of taxa being present #########
###########################################################################################
######################
def est_prob(source,induced,whole,N):
	'''Use distribution factorization from random markov field to estimate the probability of
	seeing the subset of taxa in the sample. Takes in the ROW NUMBERS of the represented edges,
	 along with the whole network.'''
	##use the whole network to get a normalizing constant...on the other hand for many purposes
	#this doesn't matter....
	num_tot_edges = int(len(whole))
	all_freq = array([whole['source_freq'][j] for j in range(0,num_tot_edges)])
	s_freq = delete(all_freq,range(1,num_tot_edges,2))
	t_freq = delete(all_freq,range(0,num_tot_edges,2))
	all_st_freq = array([whole['weight'][j] for j in range(0,num_tot_edges)])
	all_p_s_given_t = all_st_freq/all_freq
	p_s_given_t = delete(all_p_s_given_t,range(1,num_tot_edges,2))
	p_t_given_s = delete(all_p_s_given_t,range(0,num_tot_edges,2))
	normalizer = ones(num_tot_edges//2) + s_freq/N + t_freq/N + 0.5*p_s_given_t + 0.5*p_t_given_s
	Z = prod(normalizer)
	#now use the ROW NUMBERS OF THE represented edges and present edges
	edge_there = argwhere([node in induced[0] for node in source[0]])
	repped = delete(source,edge_there)
	uni_induced = delete(induced,range(1,len(induced[0]),2))
# 	print(len(uni_induced))
# 	print(len(repped))                 
	psi1 = log([0.5*all_p_s_given_t[i] + 0.5*all_p_s_given_t[i+1] for i in uni_induced])
	psi2 = log([all_freq[i]/N for i in repped])
	produ = sum(psi1) + sum(psi2) - log(Z)
	return produ

def find_cliques(node, network):
	'''find all the maximal cliques containing a given node, where network is an adjacency graph
	given as a numpy array'''
	adja = nonzero(network[node])[0]
	lvls = [[[i,node] for i in adja]]
	lvl_num = 0
	children = len(adja)
	comp_cliqs = []
	#create a tree so that tracing up from each leaf gives a maximal clique including the 
	#node in question, and grow cliques
	while children > 0:
		this_lvl = []
		this_lvl_kids = []
		for child in lvls[lvl_num]:
			lineage = child[1:]
			siblings_ind = where([chld[1:] == lineage for chld in lvls[lvl_num]])[0]
			siblings = [lvls[lvl_num][sind][0] for sind in siblings_ind]
			offspring_ind = nonzero(network[child[0],siblings])[0]
			offspring = [siblings[off_in] for off_in in offspring_ind] 
			this_lvl_kids = this_lvl_kids + [len(offspring)]
			this_lvl = this_lvl + [[off] + child for off in offspring]
			if len(offspring) == 0:
				comp_cliqs = comp_cliqs + [child]
		children = max(this_lvl_kids)
		lvl_num = lvl_num+1
		lvls = lvls + [this_lvl]
	cl_sizes = [len(cl) for cl in comp_cliqs]
	cl_sz = unique(cl_sizes)
	cliques = []
	for sz in cl_sz:
		sz_cliques = where(cl_sizes == sz)[0]
		cliques = cliques + [unique([sort(comp_cliqs[ind]) for ind in sz_cliques])]
	return cliques

def psi_over_psi(k1,k2,cliq, rm = 0.8):
	'''calculate ratio between psi(k+1) & psi(k)'''
	diff = abs(k1 - k2)
	k = min(k1,k2)
	n = len(cliq)
	fis = zeros(diff)
	rg = range(diff)
	for i in rg:
		fis[i] = (2*(1-rm)/(n-1))*(k+i) + rm
	f = prod(fis)
	if k1>k2:
		return f
	else:
		return 1/f

def diff_cliques(s1,s2,network):
	''''Take two different samples and identify which cliques are different. Then, compute the
	ratio of probabilities between the two configurations based on the rule that 
	psi(k+1)/psi(k) depends on whether or not k is above or below half the size of the clique.'''
	diff = s1 - s2
	diff_nodes = nonzero(diff)[0]
	rel_cliques = [find_cliques(node, network) for node in diff_nodes]
	rel_cliques = [i for sl in rel_cliques for i in sl]
	cl_len = [len(cl) for cl in rel_cliques]
	cl_sz = unique(cl_len)
	the_rel_cliques = []
	for sz in cl_sz:
		sz_cliques = where(cl_len == sz)[0]
		this_sz = [rel_cliques[ibd] for ibd in sz_cliques]
		this_sz = unique(this_sz,axis = 0)
		the_rel_cliques = the_rel_cliques + [i for i in this_sz]
	s1pres = nonzero(s1)[0]
	s2pres = nonzero(s2)[0]
	k1s = [sum([int(taxa in cliq) for taxa in s1pres]) for cliq in the_rel_cliques]
	k2s = [sum([int(taxa in cliq) for taxa in s2pres]) for cliq in the_rel_cliques]
	comparison = [psi_over_psi(k1s[i],k2s[i],the_rel_cliques[i]) for i in range(len(the_rel_cliques))]
	comparison = prod(comparison)
	return comparison
	
def diffusion_ivp(known_on,known_off, network, suspected =0.5, non_suspected = 0, probably = 1, sample = [], all = False):
	'''Using diffusion on the graph, rank the nodes we don't know about. Solve diffusion with
	initial condition being 1 for known on nodes and -1 for known off nodes. Network should be
	the adjacency graph (probably unweighted) given as a numpy array.'''
	#construct graph laplacian
	L = diag(sum(network,axis = 0)) - network
	#then its easy
	[eval, evec] = eig(L)
	if len(sample) == 0:
		#set up initial conditions, giving some initial mass to the nodes we are questioning about
		u0 = suspected*ones(len(network))
		u0[known_on] = probably #option to give less weight to ones we thing are on
		u0[known_off] = non_suspected #option to give some weight to the ones we think are off.
	elif len(sample) == len(network):
		u0 = sample.astype(float)
	else:
		print('Initial distribution must cover whole network')
		return False
	u0 = u0/norm(u0)
	#find the coefficients c_i
	ci = solve(evec,u0)
	zers = where(abs(eval) > 10**(-15)) #find nonzero eigenvalues
	eqs = delete(range(len(eval)),zers) #and zero eigenvalues
	equib = dot(evec[:,eqs],ci[eqs]) #equilibrium solution to the diffusion(the kernel of L)
	rel_c = ci[zers]
	rel_l = eval[zers]
	ratio_c_l = rel_c/rel_l #compute the transient behavior 
	known = known_on + known_off
	known = array(known)
	if all:
		unknown = array(range(len(network)))
	else:
		unknown = delete(array(range(len(network))),known)#nodes that we don't know about (on/off)
	equib_unk = equib[unknown].round(6)
	strengths = array([dot(ratio_c_l,transpose(evec[i,zers])) for i in unknown]).flatten()
	strengths = real(strengths.round(6))
	tmp = equib_unk.argsort()
	ranked = empty(len(tmp),int)
	requib = empty(len(tmp))
	rstren = empty(len(tmp))
	requib = equib_unk[tmp[::-1]]
	rstren = list(strengths[tmp[::-1]])
	ranked = array(unknown)[tmp[::-1]]
	ranked = list(ranked)
	srt_str = strengths[tmp[::-1]]
	ust, whch = unique(equib_unk, return_inverse = True)
	if len(ust) != len(equib_unk):
		for j in range(len(ust)):
			if sum([j == w for w in whch]) > 1:
				tied_guys = unknown[where(equib_unk == ust[j])]
				tied_which = [ranked.index(guy) for guy in tied_guys]
				their_str = srt_str[tied_which]
				tmp2 = their_str.argsort()
				them_rked = empty(len(tied_guys))
				them_rked = list(tied_guys[tmp2[::-1]])
				their_str_rk = empty(len(tmp2))
				their_str_rk = their_str[tmp2[::-1]]
				ust2,wch2 = unique(their_str,return_inverse = True)
				if len(ust2) != len(their_str):
					for k in range(len(ust2)):
						if sum([k==ww for ww in wch2])>1:
							stl_tied = tied_guys[where(their_str == ust2[k])]
							stl_which = [them_rked.index(guyy) for guyy in stl_tied]
							stl_rank = min(stl_which)
							for stl in stl_tied:
								them_rked.remove(stl)
							them_rked.insert(stl_rank,list(stl_tied))
				tied_rank = min(tied_which)
				for po in them_rked:	
					if shape(po) != ():
						for ppo in po:
							ranked.remove(ppo)
					else:
						ranked.remove(po)
				ranked.insert(tied_rank,list(them_rked))
				for po2 in their_str_rk:
					rstren.remove(po2)
				rstren.insert(tied_rank,list(their_str_rk))
		return [ranked, real(requib), real(rstren)]
	else:
		return [ranked,real(requib), real(rstren)]
	
def diffusion_bdvp(known_on,known_off, network):
	'''Using diffusion on the graph, rank the nodes we don't know about. Find equilibrium
	of solution with "boundary values" given by known node values (0,1). Network should be
	the adjacency graph (probably unweighted) given as a numpy array.'''
	#construct graph laplacian
	L = diag(sum(network,axis = 0)) - network
	known = known_on + known_off
	known = array(known)
	unknown = delete(array(range(len(network))),known)
	Lrestr = L[unknown]
	Lrestr = Lrestr[:,unknown]
	force = -array([sum(L[i,known_on]) for i in unknown])
	equib = array(lstsq(Lrestr,force))[0]
	tmp = equib.argsort()
	ranked = empty(len(tmp),int)
	requib = empty(len(tmp),int)
	requib = equib[tmp[::-1]]
	ranked = array(unknown)[tmp[::-1]]
	ranked = list(ranked)
	ust, whch = unique(equib.round(6), return_inverse = True)
	if len(ust) != len(equib):
		for j in range(len(ust)):
			if sum([j == w for w in whch]) > 1:
				tied_guys = unknown[where(equib.round(6) == ust[j])]
				tied_which = where([guy in tied_guys for guy in ranked])
				tied_rank = min(tied_which[0])
				ranked = list(ranked)
				for gy in tied_guys:
					ranked.remove(gy)
				ranked.insert(tied_rank,list(tied_guys))
		return [ranked,real(requib)]
	else:
		return [ranked,real(requib)]

	
def diffusion_forced(known_on,known_off, network):
	'''Using diffusion on the graph, rank the nodes we don't know about. Find equilibrium
	of solution with forcing (pm 1) on the known nodes. Network should be
	the adjacency graph (probably unweighted) given as a numpy array.'''
	#construct graph laplacian
	L = diag(sum(network,axis = 0)) - network
	force = zeros(len(network))
	kon = len(known_on)
	koff = len(known_off)
	if kon == 0:
		kon = 1
	if koff == 0:
		koff = 1
	force[known_on] = 1/kon
	force[known_off] = -1/koff
	unknown = argwhere(force == 0).flatten()
	equib = array(lstsq(L,force))[0]
	rel_eq = equib[unknown]
	tmp = rel_eq.argsort()
	ranked = empty(len(tmp),int)
	requib = empty(len(tmp),int)
	requib = rel_eq[tmp[::-1]]
	ranked = array(unknown)[tmp[::-1]]
	ranked = list(ranked)
	ust, whch = unique(rel_eq, return_inverse = True)
	if len(ust) != len(rel_eq):
		for j in range(len(ust)):
			if sum([j == w for w in whch]) > 1:
				tied_guys = unknown[where(rel_eq == ust[j])]
				tied_which = where([guy in tied_guys for guy in ranked])
				tied_rank = min(tied_which[0])
				ranked = list(ranked)
				for gy in tied_guys:
					ranked.remove(gy)
				ranked.insert(tied_rank,list(tied_guys))
		return [ranked,real(requib)]
	else:
		return [ranked,real(requib)]
		
#########################################################################################

############################# Create a sample to test things with #######################
#########################################################################################

def make_sample(templ):
	'''create a random sample with the organsims in the real data'''
	rsamp = pd.DataFrame(templ['LEVEL'], columns = ['LEVEL'])
	rsamp['TAXA'] = templ['TAXA']
	N = len(rsamp)
	abds = zeros(N)
	min20N = min([N,63])
	nonzer = sample(range(N), k = min20N)
	nonzerval = rand(len(nonzer))
	abds[nonzer] = nonzerval
	rsamp['abundance'] = abds
	return rsamp

def get_sample(daata):
	'''grab a real sample (column) from the template data'''
	samp = pd.DataFrame(daata['LEVEL'], columns = ['LEVEL'])
	samp['TAXA'] = daata['TAXA']
	L = len(daata.columns) - 2
	this = randint(L)+2
	which = daata.columns[this]
	samp['abundance'] = daata[which]
	samp_stype = daata.columns[this]
	return [samp,samp_stype]
	
#########################################################################################


########Flat it
def flat_two_deep(li):
	firstpass = []
	for item in li:
		if isinstance(item,list):
			firstpass = firstpass + item
		else:
			firstpass = firstpass + [item]
	final = []
	for item2 in firstpass:
		if isinstance(item2,list):
			final = final + item2
		else:
			final = final + [item2]
	return final
			







	