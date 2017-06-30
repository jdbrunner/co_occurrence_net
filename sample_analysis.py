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
from co_occ_funs import *


######################
sample_name = sys.argv[1]
# coin_net_name = sys.argv[2]
# coin_mat_name = sys.argv[3]
# coocc_net_name = sys.argv[4]
# coocc_mat_name = sys.argv[5]
# level = sys.argv[6]

coocc_net_name = sys.argv[2]
coocc_mat_name = sys.argv[3]
level = sys.argv[4]
sample_name_short = sys.argv[5]

sample = pd.read_csv(sample_name, sep = ' ')
# coin_net = pd.read_csv(coin_net_name, sep = '\t')
# coin_mat = pd.read_csv(coin_mat_name, sep = '\t', index_col = 0)
coocc_net = pd.read_csv(coocc_net_name, sep = '\t')
coocc_mat = pd.read_csv(coocc_mat_name, sep = '\t', index_col = 0)

#trim the sample to just the level we are looking at
sample = sample.iloc[where([i == level for i in sample['LEVEL']])]

wrong_guy = 'Paracoccus' 


#who is there
N = len(sample)
present = nonzero(sample['abundance'])[0]
present_names = sample['TAXA'][present]


#and trim off Taxa that aren't in my network
in_net = [(ta in coocc_mat.index) for ta in sample['TAXA']]
#but first get a list of Taxa not in my network that are here, just in case
not_in = sample.iloc[where(invert(in_net))]
pres_not_in = nonzero(not_in['abundance'])[0]
pres_not_in_names = not_in['TAXA'][pres_not_in]
#trim
sample_trim = pd.DataFrame(sample.iloc[where(in_net)].values, columns = sample.columns)


net = coocc_mat.values
mean_samp = mean(sample_trim['abundance'])
samp_on = list(where(sample_trim['abundance'].values > mean_samp)[0])
samp_off = list(where(sample_trim['abundance'].values == 0)[0])

ivp_r = empty(len(sample_trim))
bdvp_r = empty(len(sample_trim))
forced_r = empty(len(sample_trim))


ivp = diffusion_ivp(samp_on,samp_off,net,sample = sample_trim['abundance'].values, all = True)
bdvp = diffusion_bdvp(samp_on,samp_off,net)
forced = diffusion_forced(samp_on,samp_off,net)

print(ivp[0])
print(bdvp[0])
print(forced[0])



rank_rg_ivp = len(ivp[0]) 
for li in ivp[0]:
	if isinstance(li,list):
		rank_rg_ivp = rank_rg_ivp+len(li) - 1

rank_rg_bdvp = len(bdvp[0]) 
rank_rg_forced = len(forced[0]) 

ivp_cols = color_picker2(range(-1,rank_rg_ivp + 1), the_map = cm.coolwarm)[::-1]
bdvp_cols = color_picker2(range(-1,rank_rg_bdvp + 1), the_map = cm.coolwarm)[::-1]
forced_cols = color_picker2(range(-1,rank_rg_forced + 1), the_map = cm.coolwarm)[::-1]

ivp_cold = array([' '*len(ivp_cols[0])]*len(sample_trim))
bdvp_cold = array([' '*len(ivp_cols[0])]*len(sample_trim))
forced_cold = array([' '*len(ivp_cols[0])]*len(sample_trim))

ivp_r[samp_off] = rank_rg_ivp
bdvp_r[samp_off] = rank_rg_bdvp
forced_r[samp_off] = rank_rg_forced

ivp_cold[samp_off] = ivp_cols[-1]
bdvp_cold[samp_off] = bdvp_cols[-1]
forced_cold[samp_off] = forced_cols[-1]

ivp_r[samp_on] = -1
bdvp_r[samp_on] = -1
forced_r[samp_on] = -1

ivp_cold[samp_on] = ivp_cols[0]
bdvp_cold[samp_on] = bdvp_cols[0]
forced_cold[samp_on] = forced_cols[0]

l1 = 0
for rk1 in ivp[0]:
	if isinstance(rk1, list):
		for rrk1 in rk1:
			ivp_r[rrk1] = l1
			ivp_cold[rrk1] = ivp_cols[l1+1]
			l1 = l1+1
	else:
		ivp_r[rk1] = l1
		ivp_cold[rk1] = ivp_cols[l1+1]
		l1 = l1+1
	

l2 = 0
for rk2 in bdvp[0]:
	bdvp_r[rk2] = l2
	bdvp_cold[rk2] = ivp_cols[l2+1]
	l2 = l2+1
	
l3 = 0
for rk3 in forced[0]:
	forced_r[rk3] = l3
	forced_cold[rk3] = forced_cols[l3+1]
	l3 = l3+1

sample_trim['ivp'] = ivp_r
sample_trim['bdvp'] = bdvp_r
sample_trim['forced'] = forced_r

sample_trim['ivp_col'] = ivp_cold
sample_trim['bdvp_col'] = bdvp_cold
sample_trim['forced_col'] = forced_cold


ivp_r_net = -ones(len(coocc_net))
bdvp_r_net = -ones(len(coocc_net))
forced_r_net = -ones(len(coocc_net))

ivp_col_net = array(len(coocc_net)*[ivp_cols[0]])
bdvp_col_net = array(len(coocc_net)*[bdvp_cols[0]])
forced_col_net = array(len(coocc_net)*[forced_cols[0]])

 
order_in_net = [where(coocc_net['source'] == tax)[0] for tax in sample_trim['TAXA']]
q1 = 0
for lo in order_in_net:
	ivp_r_net[lo] = sample_trim['ivp'].values[q1]
	ivp_col_net[lo] = sample_trim['ivp_col'].values[q1]
	bdvp_r_net[lo] = sample_trim['bdvp'].values[q1]
	bdvp_col_net[lo] = sample_trim['bdvp_col'].values[q1]
	forced_r_net[lo] = sample_trim['forced'].values[q1]
	forced_col_net[lo] = sample_trim['forced_col'].values[q1]
	q1 = q1 + 1
	

coocc_net[sample_name_short+'_ivp'] = ivp_r_net
coocc_net[sample_name_short+'_bdvp'] = bdvp_r_net
coocc_net[sample_name_short+'_forced'] = forced_r_net

coocc_net[sample_name_short+'_ivp_col'] = ivp_col_net
coocc_net[sample_name_short+'_bdvp_col'] = bdvp_col_net
coocc_net[sample_name_short+'_forced_col'] = forced_col_net

coocc_net.to_csv(coocc_net_name, sep = '\t')

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


 