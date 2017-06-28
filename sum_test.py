from pylab import *
import pandas as pd

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
	[eval, evec] = eig(-L)
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
	requib = empty(len(tmp),int)
	requib = rel_eq[tmp[::-1]]
	ranked = array(unknown)[tmp[::-1]]
	ust, whch = unique(rel_eq, return_inverse = True)
	if len(ust) != len(rel_eq):
		ties = []
		print(ust)
		for j in range(len(ust)):
			if sum([j == w for w in whch]) > 1:
				ties = ties + [unknown[where(rel_eq == ust[j])]]
		return [ranked,ties,requib]
	else:
		return [ranked,requib]
	
	
samp1 = zeros(7)
samp2 = zeros(7)
samp1[[0,1,2,3]] = 1
samp2[[0,1,2,4]] = 1	

	
adj = zeros([7,7])
adj[3,2] = 1
adj[5,6] = 1
adj[1,2] = 1
adj[0,3] = 1
adj[0,1] = 1
adj[4,2] = 1
adj[3,1] = 1
adj[0,2] = 1
adj[1,4] = 1
adj[3,5] = 1
adj[4,6] = 1
adj = adj + transpose(adj)



# print(diff_cliques(samp1,samp2,adj))
print(diffusion_ivp([0,1,2],[6],adj))
print(diffusion_bdvp([0,1,2],[6],adj))
print(diffusion_forced([0,1,2],[6],adj))
# 
adj2 = array([[0,1,1,0],[1,0,1,0],[1,1,0,1],[0,0,1,0]])
# 
print(diffusion_ivp([0,2],[],adj2))
print(diffusion_bdvp([0,2],[],adj2))
print(diffusion_forced([0,2],[],adj2))

print(diffusion_ivp([2],[0],adj2))
print(diffusion_bdvp([2],[0],adj2))
print(diffusion_forced([2],[0],adj2))















