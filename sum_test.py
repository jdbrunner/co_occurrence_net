from pylab import *
import pandas as pd
from co_occ_funs import *

	
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



#print(diff_cliques(samp1,samp2,adj))
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















