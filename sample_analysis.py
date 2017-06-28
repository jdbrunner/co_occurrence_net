#########################################################################################
#																						#		
#				Sample Analysis															#
#																						#
#										Jim Brunner 									#
#																						#
#																						#
#																						#
#########################################################################################
#
# We want to use the co-occurence network to say something about a
#	sample
#
#
## modules needed
from pylab import *
import pandas as pd
import sys
import matplotlib.pyplot as plt

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
	print(len(uni_induced))
	print(len(repped))                 
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
	
def diffusion_ivp(known_on,known_off, network):
	'''Using diffusion on the graph, rank the nodes we don't know about. Solve diffusion with
	initial condition being 1 for known on nodes and -1 for known off nodes. Network should be
	the adjacency graph (probably unweighted) given as a numpy array.'''
	#construct graph laplacian
	L = diag(sum(network,axis = 0)) - network
	#then its easy
	[eval, evec] = eig(L)
	u0 = zeros(len(network))
	u0[known_on] = 1
	#u0[known_off] = -1
	ci = solve(evec,u0)
	zers = where(abs(eval) > 10**(-5))
	rel_c = ci[zers]
	rel_l = eval[zers]
	ratio_c_l = -rel_c/rel_l
	known = known_on + known_off
	known = array(known)
	unknown = delete(array(range(len(network))),known)
	strengths = array([dot(ratio_c_l,transpose(evec[i,zers])) for i in unknown]).flatten()
	tmp = strengths.argsort()
	ranked = empty(len(tmp),int)
	rstren = empty(len(tmp),int)
	rstren = strengths[tmp[::-1]]
	ranked = array(unknown)[tmp[::-1]]
	ust, whch = unique(strengths.round(6), return_inverse = True)
	if len(ust) != len(strengths):
		ties = []
		for j in range(len(ust)):
			if sum([j == w for w in whch]) > 1:
				ties = ties + [unknown[where(strengths.round(6) == ust[j])]]
		return [ranked,ties,rstren]
	else:
		return [ranked,rstren]
	
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
	force = array([sum(L[i,known_on]) for i in unknown])
	equib = array(solve(-Lrestr,force))
	tmp = equib.argsort()
	ranked = empty(len(tmp),int)
	requib = empty(len(tmp),int)
	requib = equib[tmp[::-1]]
	ranked = array(unknown)[tmp[::-1]]
	ust, whch = unique(equib.round(6), return_inverse = True)
	if len(ust) != len(equib):
		ties = []
		for j in range(len(ust)):
			if sum([j == w for w in whch]) > 1:
				ties = ties + [unknown[where(equib.round(6) == ust[j])]]
		return [ranked,ties,requib]
	else:
		return [ranked,requib]

	
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
	ranked = array(unknown)[tmp[::-1]]
	ust, whch = unique(rel_eq.round(6), return_inverse = True)
	if len(ust) != len(rel_eq):
		ties = []
		for j in range(len(ust)):
			if sum([j == w for w in whch]) > 1:
				ties = ties + [unknown[where(rel_eq.round(6) == ust[j])]]
		return [ranked,ties]
	else:
		return ranked
	
	


######################
sample_name = sys.argv[1]
# coin_net_name = sys.argv[2]
# coin_mat_name = sys.argv[3]
# coocc_net_name = sys.argv[4]
# coocc_mat_name = sys.argv[5]
# level = sys.argv[6]

coocc_mat_name = sys.argv[2]
level = sys.argv[3]

sample = pd.read_csv(sample_name, sep = ' ')
# coin_net = pd.read_csv(coin_net_name, sep = '\t')
# coin_mat = pd.read_csv(coin_mat_name, sep = '\t', index_col = 0)
# coocc_net = pd.read_csv(coocc_net_name, sep = '\t')
coocc_mat = pd.read_csv(coocc_mat_name, sep = '\t', index_col = 0)

#trim the sample to just the level we are looking at
sample = sample.iloc[where([i == level for i in sample['LEVEL']])]

#and trim off Taxa that aren't in my network
in_net = [(ta in coocc_mat.index) for ta in sample['TAXA']]
sample = sample.iloc[where(in_net)]

#who is there
N = len(sample)
present = nonzero(sample['abundance'])[0]
present_names = sample['TAXA'][present]
#print(len(present))

net = coocc_mat.values
mean_samp = mean(sample['abundance'])
samp_on = list(where(sample['abundance'] > mean_samp)[0])
samp_off = list(where(sample['abundance'] == 0)[0])
print(diffusion_ivp(samp_on,samp_off,net))
print(diffusion_bdvp(samp_on,samp_off,net))
#print(diffusion_forced(samp_on,samp_off,net))


#construct the network from the sample.
# samp_net = zeros([N,N])
# for ind in present:
# 	samp_net[ind] = 1
# samp_net = samp_net + transpose(samp_net)
# sample_df = pd.DataFrame(samp_net, columns = sample['TAXA'], index = sample['TAXA'])

#get the induced subgraph of the coincidence and cooccurrence networks.
# induced_coin = coin_mat[present_names]
# induced_coin = induced_coin.loc[present_names]
# induced_edges_coin = where([(coin_net['source'][i] in array(present_names)) & (coin_net['target'][i] in array(present_names))
# 			for i in range(len(coin_net))])
# induced_coin_net = coin_net.iloc[induced_edges_coin]

# induced_coocc = coocc_mat[present_names]
# induced_coocc = induced_coocc.loc[present_names]
# induced_edges_coocc = where([(coocc_net['source'][i] in array(present_names)) & (coocc_net['target'][i] in array(present_names))
# 			for i in range(len(coocc_net))])
# repped_sources = where([coocc_net['source'][i] in array(present_names) for i in range(len(coocc_net))])
# samp_prob = est_prob(repped_sources,induced_edges_coocc,coocc_net,coocc_net['num_samps'][0])
# induced_coocc_net = coocc_net.iloc[induced_edges_coocc]
# print(samp_prob)

#save the induced networks for loading into cytoscape.
# induced_coin_net.to_csv(sample_name[:-4]+'_induced_coin.tsv', sep = '\t')
# induced_coocc_net.to_csv(sample_name[:-4]+'_induced_coocc.tsv', sep = '\t')


 