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


template_name = sys.argv[1]
level = sys.argv[2]

template = pd.read_csv(template_name, sep = ' ')

def make_sample(templ):
	'''create a random sample with the organsims in the real data'''
	rsamp = pd.DataFrame(templ['LEVEL'], columns = ['LEVEL'])
	rsamp['TAXA'] = templ['TAXA']
	N = len(rsamp)
	abds = zeros(N)
	min20N = min([N,20])
	nonzer = sample(range(N), k = min20N)
	nonzerval = rand(len(nonzer))
	abds[nonzer] = nonzerval
	rsamp['abundance'] = abds
	return rsamp

rsamps = dict() 

#We want to make a network at a specific level, or set of levels
if level == 'all':
	level = array(template['LEVEL'].values)
	level = array(unique(level))
else:
	level = array([level])

#need to check that the levels given as argument actually exist
levels_allowed = array(template['LEVEL'].values)
levels_allowed = array(unique(levels_allowed))

for i in level:
	if not i in levels_allowed:
		print(levels_allowed)
		sys.exit('Bad level argument. Choose a level from the list above or type all')

for lvl in level:
		
	#start by getting the indices of the members of the level
	in_level = where(template['LEVEL'] == i)
	lvl_template = template.iloc[in_level]
	
	rsamps[lvl] = make_sample(lvl_template)
	
for i in level:
	flname1 = i+'_randsamp.txt'
	rsamps[i].to_csv(flname1, sep = ' ')


	
	
	
	
	
	
	
	
	
	
	
	