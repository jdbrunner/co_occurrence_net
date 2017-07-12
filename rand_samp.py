#########################################################################################
#																						#		
#				Random sample maker														#
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
from random import sample
from co_occ_funs import *



template_name = sys.argv[1]
level = sys.argv[2]
flder = sys.argv[3]

template = pd.read_csv(template_name, sep = ' ')


##I have to parse the level input I guess...
if level != 'all':
	lvl_list = []
	if level[0] == '[' and level[-1] == ']':
		commas = argwhere([let == ',' for let in level])
		word_start = 1
		for com in commas:
			lvl_list = lvl_list + [level[word_start:com[0]]]
			word_start = com[0]+1
		lvl_list = lvl_list+[level[word_start:-1]]
	else:
		lvl_list  = [level]
	level = lvl_list


rand_samp = True	

if rand_samp:
	rsamps = dict() 
colsamps = dict()

#We want to make a network at a specific level, or set of levels
if level == 'all':
	level = array(template['LEVEL'].values)
	level = array(unique(level))
else:
	level = array(level)

#need to check that the levels given as argument actually exist
levels_allowed = array(template['LEVEL'].values)
levels_allowed = array(unique(levels_allowed))

for i in level:
	if not i in levels_allowed:
		print(levels_allowed)
		sys.exit('Bad level argument. Choose a level from the list above or type all')

for lvl in level:
		
	#start by getting the indices of the members of the level
	in_level = where(template['LEVEL'] == lvl)
	lvl_template = template.iloc[in_level]
	if rand_samp:
		rsamps[lvl] = make_sample(lvl_template)
	
	colsamps[lvl] = get_sample(lvl_template)
	samp_len = len(colsamps[lvl][0])
	#in_samp is an array of indices (locs, not ilocs) that have non-zero abundance
	in_samp = colsamps[lvl][0].index[nonzero(colsamps[lvl][0]['abundance'])]
	not_in_samp = colsamps[lvl][0].index[where(colsamps[lvl][0]['abundance'].values == 0)]
	num_in_samp = len(in_samp)
	num_not = len(not_in_samp)
	rmv_amnt = 0#max(1,int(num_in_samp*0.05))
	rmv = randint(num_in_samp,size = rmv_amnt)
	# the indices of the taxa set to zero are in in_samp[rmv]
	colsamps[lvl][0].loc[in_samp[rmv],'abundance'] = 0
	add_amnt = 0 
	add_em = randint(num_not,size = add_amnt)
	#similar to above
	colsamps[lvl][0].loc[not_in_samp[add_em],'abundance'] = rand(add_amnt)*mean(colsamps[lvl][0]['abundance'][in_samp])
	
	colsamps[lvl][0]['Detected'] = array([False]*len(colsamps[lvl][0]))
	colsamps[lvl][0]['Removed'] = array([False]*len(colsamps[lvl][0]))
	colsamps[lvl][0]['Added'] = array([False]*len(colsamps[lvl][0]))
	
	rsamps[lvl]['Detected'] = array([False]*len(rsamps[lvl]))
	rsamps[lvl]['Removed'] = array([False]*len(rsamps[lvl]))
	rsamps[lvl]['Added'] = array([False]*len(rsamps[lvl]))
	
	det_loc_in_in_samp = delete(range(num_in_samp),rmv)
	colsamps[lvl][0].loc[in_samp[det_loc_in_in_samp],'Detected'] = True
	colsamps[lvl][0].loc[not_in_samp[add_em],'Detected'] = True
	
	colsamps[lvl][0].loc[in_samp[rmv],'Removed'] = True
	colsamps[lvl][0].loc[not_in_samp[add_em],'Added'] = True
	
	if rand_samp:
		flname1 = lvl+'/rand.tsv'
		rsamps[lvl].to_csv(flder+'/'+flname1, sep = '\t')
	
	flname2 = lvl+'/'+colsamps[lvl][1]+'.tsv'
	colsamps[lvl][0].to_csv(flder+'/'+flname2, sep = '\t')
	print(flname2)
	

	

	

	
	
	
	
	
	
	
	
	