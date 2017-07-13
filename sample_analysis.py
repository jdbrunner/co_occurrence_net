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
###pandas DataFrame usage: df['column_name'] or df.loc['row_name' (index),'column'] or df.iloc[row_number, 'column']
import sys
import matplotlib.pyplot as plt
from co_occ_funs import *


######################
#param_file = sys.argv[1]

#flder = sys.argv[5]

# with open(param_file,'r') as param:
# 	pars = list(param)[1:]

# param_starts = [par.index('=') + 2 for par in pars]
# 
# arguments = [pars[n][param_starts[n]:-1] for n in range(len(pars))]
# 
sample_name = sys.argv[1]#arguments[0]


coocc_mat_name = sys.argv[3]#arguments[1]#
node_atts_name = sys.argv[2]
level = sys.argv[4]#arguments[3]#

ntype = ''
if 'bins' in coocc_mat_name:
	ntype = '_bins'
elif 'pear' in coocc_mat_name:
	ntype = '_pear'

slsh = coocc_mat_name.rfind('/')
unders = where([letter == '_' for letter in coocc_mat_name])[0]
sample_name_short = coocc_mat_name[slsh+1:unders[-1]]+ntype #sample_name[slsh+1:-14]#sys.argv[5]#arguments[4]#

sample = pd.read_csv(sample_name, sep = '\t')
node_atts = pd.read_csv(node_atts_name, sep = '\t')
coocc_mat = pd.read_csv(coocc_mat_name, sep = '\t', index_col = 0)

#trim the sample to just the level we are looking at (shouldn't be needed with the samples I made)
samplelvl = sample.iloc[where([i == level for i in sample['LEVEL']])]

#ilocs!
rmved = where(samplelvl['Removed'])
added = where(samplelvl['Added'])

#who is there
N = len(samplelvl)
present = nonzero(samplelvl['abundance'])[0]
present_names = samplelvl['TAXA'][present]


#and trim off Taxa that aren't in my network
in_net = [(ta in coocc_mat.index) for ta in samplelvl['TAXA']]
#but first get a list of Taxa not in my network that are here, just in case
not_in = samplelvl.iloc[where(invert(in_net))]
pres_not_in = nonzero(not_in['abundance'])[0]
pres_not_in_names = not_in['TAXA'][pres_not_in]
#trim
in_net_ilocs = where(in_net)
in_net_locs = samplelvl.index[in_net_ilocs]
sample_trim = samplelvl.loc[in_net_locs]
####The above way fo defining a new DF keeps the old index!

#it's possible that anything present in the sample isn't in my network yet
### this happens also when taxa aren't connected to anything - they are cut out of the network
if shape(nonzero(sample_trim['abundance']))[1] == 0:
	sys.exit('Sample is empty or only contains unconnected taxa')

#get the adjacency matrix as an np array
net = coocc_mat.values
#average of non-zero entries
mean_samp = mean(sample_trim['abundance'].values[nonzero(sample_trim['abundance'])])
##samp_on is ilocs, not locs.
samp_on = list(where(sample_trim['abundance'].values > mean_samp)[0])
samp_off = []#list(where(sample_trim['abundance'].values == 0)[0])

#these will store the rank of the taxa in the corresponding iloc
ivp_r = empty(len(sample_trim))
bdvp_r = empty(len(sample_trim))
forced_r = empty(len(sample_trim))

##At this point, the sample ilocs need to match the adjacency matrix. This SHOULD be true, provided
#there is no bug in the trimming off of taxa not in the network. Notice, the indices will not match,
#because the indices from the full sample are kept.
ivp = diffusion_ivp(samp_on,samp_off,net,sample = sample_trim['abundance'].values, all = True)
#returns [ranked,real(requib), real(rstren)]  - if all=False, only gives for ilocs not in samp_on or samp_off
bdvp = diffusion_bdvp(samp_on,samp_off,net)
#returns [ranked,real(requib)] only gives for ilocs not in samp_on or samp_off
forced = diffusion_forced(samp_on,samp_off,net)
#returns [ranked,real(requib)] only gives for ilocs not in samp_on or samp_off


##all three give a list of ilocs ordered by rank - ties are given as lists within the list.  The first 
# gives a list of lists (tied at equilibrium), and when there are ties in the transient 
# behavior it gives 

#get the number of uniquely ranked elements in the lists
rank_rg_ivp = sum([len(eq_tie) if isinstance(eq_tie,list) else 1 for eq_tie in ivp[0]]) 
rank_rg_bdvp = len(bdvp[0]) 
rank_rg_forced = len(forced[0]) 


##Gets a list of colors to use
ivp_cols = color_picker2(range(-1,rank_rg_ivp + 1), the_map = cm.coolwarm)[::-1]
bdvp_cols = color_picker2(range(-1,rank_rg_bdvp + 1), the_map = cm.coolwarm)[::-1]
forced_cols = color_picker2(range(-1,rank_rg_forced + 1), the_map = cm.coolwarm)[::-1]
samp_cols = color_picker2(range(len(unique(sample_trim['abundance'].values))), the_map = cm.coolwarm)

#the transient behavior flattened.
flat_transient = real(array([val for sublist in ivp[2] for val in sublist]))

##Gets a list of colors to use
ivp_eq_cols = color_picker2(ivp[1],the_map = cm.coolwarm, weighted = True)
ivp_tr_cols = color_picker2(flat_transient,the_map = cm.coolwarm, weighted = True)
bdvp_eq_cols = color_picker2(bdvp[1],the_map = cm.coolwarm, weighted = True)
forced_eq_cols = color_picker2(forced[1],the_map = cm.coolwarm, weighted = True)

### need slots for hex values
ivp_cold = array([' '*len(ivp_cols[0])]*len(sample_trim))
bdvp_cold = array([' '*len(ivp_cols[0])]*len(sample_trim))
forced_cold = array([' '*len(ivp_cols[0])]*len(sample_trim))

#going to store the actual values (not just the list ordered by rank)
ivp_vals_eq = empty(len(sample_trim))
ivp_vals_tr = empty(len(sample_trim))
bdvp_vals = empty(len(sample_trim))
forced_vals = empty(len(sample_trim))

### need slots for hex values
ivp_vals_eq_col = array([' '*len(ivp_cols[0])]*len(sample_trim))
ivp_vals_tr_col = array([' '*len(ivp_cols[0])]*len(sample_trim))
bdvp_vals_col = array([' '*len(ivp_cols[0])]*len(sample_trim))
forced_vals_col = array([' '*len(ivp_cols[0])]*len(sample_trim))

##the off nodes automatically get the highest rank (they are automatically last)
ivp_r[samp_off] = rank_rg_ivp
bdvp_r[samp_off] = rank_rg_bdvp
forced_r[samp_off] = rank_rg_forced

#and the off nodes are 0 by def for all by the forcing.
ivp_vals_eq[samp_off] = 0
ivp_vals_tr[samp_off] = 0
bdvp_vals[samp_off] = 0
forced_vals[samp_off] = min(forced[1])

#off gets the coolest color
ivp_vals_eq_col[samp_off] = ivp_eq_cols[0]
ivp_vals_tr_col[samp_off] = ivp_tr_cols[0]
bdvp_vals_col[samp_off] = bdvp_eq_cols[0]
forced_vals_col[samp_off] = forced_eq_cols[0]


ivp_cold[samp_off] = ivp_cols[-1]
bdvp_cold[samp_off] = bdvp_cols[-1]
forced_cold[samp_off] = forced_cols[-1]

#on nodes are first automatically
ivp_r[samp_on] = -1
bdvp_r[samp_on] = -1
forced_r[samp_on] = -1

ivp_vals_eq[samp_on] = max(ivp[1])
ivp_vals_tr[samp_on] = max(flat_transient)
bdvp_vals[samp_on] = max(bdvp[1])
forced_vals[samp_on] = max(forced[1])

ivp_vals_eq_col[samp_on] = ivp_eq_cols[-1]
ivp_vals_tr_col[samp_on] = ivp_tr_cols[-1]
bdvp_vals_col[samp_on] = bdvp_eq_cols[-1]
forced_vals_col[samp_on] = forced_eq_cols[-1]


ivp_cold[samp_on] = ivp_cols[0]
bdvp_cold[samp_on] = bdvp_cols[0]
forced_cold[samp_on] = forced_cols[0]

#ivp returns a list of lists (of lists if there are ties down to the transient behavior)
l1 = 0
l11 = 0
p1 = 0
for rk1 in ivp[0]:
	if isinstance(rk1, list):
		for rrk1 in rk1:
			#ties at the transient level should get the same rank
			ivp_r[rrk1] = l1
			ivp_cold[rrk1] = ivp_cols[l11+1]
			#equib values is a flat list
			ivp_vals_eq[rrk1] = ivp[1][p1]
			ivp_vals_eq_col[rrk1] = ivp_eq_cols[p1]
			if isinstance(rrk1, list):
				le = len(rrk1)
				cols = ivp_tr_cols[p1:p1+le]
			else:
				le = 1
				cols = ivp_tr_cols[p1]
			l1 = l1+le
			ivp_vals_tr[rrk1] = flat_transient[p1:p1+le]
			ivp_vals_tr_col[rrk1] = cols
			p1 = p1 + le
			l11 = l11 + 1
	else:
		ivp_r[rk1] = l1
		ivp_cold[rk1] = ivp_cols[l11+1]
		ivp_vals_eq[rrk1] = ivp[1][p1]
		ivp_vals_eq_col[rk1] = ivp_eq_cols[p1]
		if isinstance(rk1, list):
			le = len(rk1)
			ivp_tr_cols[p1:p1+le]
		else:
			le = 1
			ivp_tr_cols[p1]
		l1 = l1+le
		l11 = l11 +1
		ivp_vals_tr[rk1] = flat_transient[p1:p1+le]
		ivp_vals_tr_col[rk1] = cols
		p1 = p1 + le
	
#bdvp has a list of lists (if there are ties): [t1,t2,[t3,t4],t5] means t3 & t4 are tied
l2 = 0
l22 = 0
p2 = 0
for rk2 in bdvp[0]:
	bdvp_r[rk2] = l2
	bdvp_cold[rk2] = bdvp_cols[l22+1]
	bdvp_vals[rk2] = bdvp[1][p2]
	bdvp_vals_col[rk2] = bdvp_eq_cols[p2]
	if isinstance(rk2, list):
		le = len(rk2)
	else:
		le = 1
	p2 = p2 + le
	l2 = l2+le
	l22 = l22 + 1

	
#again, a list of lists (if there are ties)
l3 = 0
l33 = 0
p3 = 0
for rk3 in forced[0]:
	forced_r[rk3] = l3
	forced_cold[rk3] = forced_cols[l33+1]
	forced_vals[rk3] = forced[1][p3]
	forced_vals_col[rk3] = forced_eq_cols[p3]
	if isinstance(rk3, list):
		le = len(rk3)
	else:
		le = 1
	p3 = p3+le
	l3 = l3+le
	l33 = l33 + 1

sample_trim['ivp'+'_'+sample_name_short] = ivp_r
sample_trim['bdvp'+'_'+sample_name_short] = bdvp_r
sample_trim['forced'+'_'+sample_name_short] = forced_r
sample_trim['ivp_eq'+'_'+sample_name_short] = ivp_vals_eq
sample_trim['ivp_transient'+'_'+sample_name_short] = ivp_vals_tr
sample_trim['bdvp_eq'+'_'+sample_name_short] = bdvp_vals
sample_trim['forced_eq'+'_'+sample_name_short] = forced_vals

sample_trim['ivp_col'+'_'+sample_name_short] = ivp_cold
sample_trim['bdvp_col'+'_'+sample_name_short] = bdvp_cold
sample_trim['forced_col'+'_'+sample_name_short] = forced_cold
samp_ord = unique(sample_trim['abundance'].values,return_inverse = True)[1]
sample_trim['sample_col'] = array(samp_cols)[samp_ord.astype(int)]
sample_trim['ivp_eq_color'+'_'+sample_name_short] = ivp_vals_eq_col
sample_trim['ivp_transient_color'+'_'+sample_name_short] = ivp_vals_tr_col
sample_trim['bdvp_eq_color'+'_'+sample_name_short] = bdvp_vals_col
sample_trim['forced_eq_color'+'_'+sample_name_short] = forced_vals_col


rmvd_ilocs = where([rtax in sample.loc[sample.index[rmved],'TAXA'].values for rtax in sample_trim['TAXA']])
rmvd_locs = sample_trim.index[rmvd_ilocs]

###This will print the ranks of the taxa removed from the sample (to test) to the terminal
term_out = False
if term_out:
	for rmd in rmvd_locs:
		print(sample_trim.loc[rmd,'ivp'+'_'+sample_name_short])
		num_beat = sum((sample_trim['ivp'+'_'+sample_name_short].values > sample_trim.loc[rmd,'ivp']).astype(int))
		print(num_beat)

	for rmd in rmvd_locs:
		print(sample_trim.loc[rmd,'bdvp'+'_'+sample_name_short])
		num_beat = sum((sample_trim['bdvp'+'_'+sample_name_short].values > sample_trim.loc[rmd,'bdvp']).astype(int))
		print(num_beat)


	for rmd in rmvd_locs:
		print(sample_trim.loc[rmd,'forced'+'_'+sample_name_short])
		num_beat = sum((sample_trim['forced'+'_'+sample_name_short].values > sample_trim.loc[rmd,'forced']).astype(int))
		print(num_beat)


##Save the sample with new columns. And clean it up a bit, make it nice for cytoscape

max_r = len(sample_trim)

sample['ivp'+'_'+sample_name_short] = max_r*ones(len(sample))
sample['ivp_eq'+'_'+sample_name_short] = zeros(len(sample))
sample['ivp_transient'+'_'+sample_name_short] = zeros(len(sample))
sample['sample_col_'+sample_name_short] = array([' ']*len(sample))

sample.index = sample['TAXA']
sample.loc[sample_trim['TAXA'],'ivp'+'_'+sample_name_short] = sample_trim['ivp'+'_'+sample_name_short].values
sample.loc[sample_trim['TAXA'],'ivp_eq'+'_'+sample_name_short] = sample_trim['ivp_eq'+'_'+sample_name_short].values
sample.loc[sample_trim['TAXA'],'ivp_transient'+'_'+sample_name_short] = sample_trim['ivp_transient'+'_'+sample_name_short].values
sample.loc[sample_trim['TAXA'],'sample_col_'+sample_name_short] = sample_trim['sample_col'].values


other_probs = False
if other_probs:
	sample_trim['bdvp'+'_'+sample_name_short] = bdvp_r
	sample_trim['forced'+'_'+sample_name_short] = forced_r
	sample_trim['bdvp_eq'+'_'+sample_name_short] = bdvp_vals
	sample_trim['forced_eq'+'_'+sample_name_short] = forced_vals
	
	sample['bdvp'+'_'+sample_name_short] = max_r*ones(len(sample))
	sample['forced'+'_'+sample_name_short] = max_r*ones(len(sample))
	sample['bdvp_eq'+'_'+sample_name_short] = zeros(len(sample))
	sample['forced_eq'+'_'+sample_name_short] = zeros(len(sample))
	
	sample.loc[sample_trim['TAXA'],'bdvp'+'_'+sample_name_short] = sample_trim['bdvp'+'_'+sample_name_short]	
	sample.loc[sample_trim['TAXA'],'forced'+'_'+sample_name_short] = sample_trim['forced'+'_'+sample_name_short]
	sample.loc[sample_trim['TAXA'],'bdvp_eq'+'_'+sample_name_short] = sample_trim['bdvp_eq'+'_'+sample_name_short]
	sample.loc[sample_trim['TAXA'],'forced_eq'+'_'+sample_name_short] = sample_trim['forced_eq'+'_'+sample_name_short]

color_it = True
if color_it:
	sample_trim['ivp_col'+'_'+sample_name_short] = ivp_cold
	
	samp_ord = unique(sample_trim['abundance'].values,return_inverse = True)[1]
	sample_trim['sample_col'] = array(samp_cols)[samp_ord.astype(int)]
	sample_trim['ivp_eq_color'+'_'+sample_name_short] = ivp_vals_eq_col
	sample_trim['ivp_transient_color'+'_'+sample_name_short] = ivp_vals_tr_col
	
	sample['ivp_col'+'_'+sample_name_short] =  array([' ']*len(sample))
	sample['ivp_eq_color'+'_'+sample_name_short] =  array([' ']*len(sample))
	sample['ivp_transient_color'+'_'+sample_name_short] =  array([' ']*len(sample))
	
	sample.loc[sample_trim['TAXA'],'ivp_col'+'_'+sample_name_short] = sample_trim['ivp_col'+'_'+sample_name_short].values
	sample.loc[sample_trim['TAXA'],'ivp_eq_color'+'_'+sample_name_short] = sample_trim['ivp_eq_color'+'_'+sample_name_short].values
	sample.loc[sample_trim['TAXA'],'ivp_transient_color'+'_'+sample_name_short] = sample_trim['ivp_transient_color'+'_'+sample_name_short].values

	if other_probs:
		sample_trim['bdvp_col'+'_'+sample_name_short] = bdvp_cold
		sample_trim['forced_col'+'_'+sample_name_short] = forced_cold
		sample_trim['bdvp_eq_color'+'_'+sample_name_short] = bdvp_vals_col
		sample_trim['forced_eq_color'+'_'+sample_name_short] = forced_vals_col
	
		sample['bdvp_col'+'_'+sample_name_short] =  array([' ']*len(sample))
		sample['forced_col'+'_'+sample_name_short] =  array([' ']*len(sample))
		sample['bdvp_eq_color'+'_'+sample_name_short] =  array([' ']*len(sample))
		sample['forced_eq_color'+'_'+sample_name_short] =  array([' ']*len(sample))
	
		sample.loc[sample_trim['TAXA'],'bdvp_col'+'_'+sample_name_short] = sample_trim['bdvp_col'+'_'+sample_name_short]
		sample.loc[sample_trim['TAXA'],'forced_col'+'_'+sample_name_short] = sample_trim['forced_col'+'_'+sample_name_short]
		sample.loc[sample_trim['TAXA'],'bdvp_eq_color'+'_'+sample_name_short] = sample_trim['bdvp_eq_color'+'_'+sample_name_short]
		sample.loc[sample_trim['TAXA'],'forced_eq_color'+'_'+sample_name_short] = sample_trim['forced_eq_color'+'_'+sample_name_short]

##Let's add some node attributes to the sample file
node_atts.index = node_atts['source']
sample['community_cluster_'+sample_name_short] = node_atts.loc[sample['TAXA'],'commun_cluster']
sample['spectral_cluster_'+sample_name_short] = node_atts.loc[sample['TAXA'],'spect_cluster']
sample['common_location_'+sample_name_short] = node_atts.loc[sample['TAXA'], 'source_type']



#sample_trim.drop('LEVEL',axis = 1, inplace = True)
for col in unique(sample.columns.values):
	if 'Unnamed' in col:
		sample.drop(col,axis = 1, inplace = True)
	if 'TAXA.1' in col:
		sample.drop(col,axis = 1, inplace = True)

sample.to_csv(sample_name, sep = '\t')




 