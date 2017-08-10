import pandas as pd
from co_occ_funs import *
import os
from pylab import *

merdat = pd.read_csv('merged2.txt', sep = ' ')
spec_lv = merdat[merdat.LEVEL == 'species']
gen_lv = merdat[merdat.LEVEL == 'genus']

# gender = pd.read_csv('gender.tsv',sep = '\t', index_col = 0)
# 
# gender.Gender[[isinstance(ent,float) for ent in gender.Gender]] = 'missing'

ids = merdat.columns[2:].values
mal_fem = array(['Female' if 'female' in lab else 'Male' for lab in ids])
mal_fem[where(invert(['male' in lab for lab in ids]))] = 'missing'

gender = pd.DataFrame(mal_fem, columns = ['Gender'], index = ids)
gender['Name'] = ['nm']*len(gender)


gen_pears = array(os.listdir('08_08_17_09_networks/genus/pears'))
gen_bins = array(os.listdir('08_08_17_09_networks/genus/bins'))
spec_pears = array(os.listdir('08_08_17_09_networks/species/pears'))
spec_bins = array(os.listdir('08_08_17_09_networks/species/bins'))


gen_p_edges = gen_pears[where(['list' in flnm for flnm in gen_pears])]
gen_p_ndat = gen_pears[where(['node_data' in flnm for flnm in gen_pears])]

gen_b_edges = gen_bins[where(['list' in flnm for flnm in gen_bins])]
gen_b_ndat = gen_bins[where(['node_data' in flnm for flnm in gen_bins])]

spec_p_edges = spec_pears[where(['list' in flnm for flnm in spec_pears])]
spec_p_ndat = spec_pears[where(['node_data' in flnm for flnm in spec_pears])]

spec_b_edges = spec_bins[where(['list' in flnm for flnm in spec_bins])]
spec_b_ndat = spec_bins[where(['node_data' in flnm for flnm in spec_bins])]

for net in range(len(gen_p_edges)):
	edges = pd.read_csv('08_08_17_09_networks/genus/pears/'+gen_p_edges[net], sep = '\t')
	node_dat = pd.read_csv('08_08_17_09_networks/genus/pears/'+gen_p_ndat[net], sep = '\t')
	with_gen = make_meta_from_file(edges, gender, gen_lv, existing_table = node_dat)
	if len(with_gen)> 0:
		with_gen.to_csv('08_08_17_09_networks/genus/pears/'+gen_p_ndat[net], sep = '\t', index = False)

for net in range(len(gen_b_edges)):
	edges = pd.read_csv('08_08_17_09_networks/genus/bins/'+gen_b_edges[net], sep = '\t')
	node_dat = pd.read_csv('08_08_17_09_networks/genus/bins/'+gen_b_ndat[net], sep = '\t')
	with_gen = make_meta_from_file(edges, gender, gen_lv, existing_table = node_dat)
	if len(with_gen)> 0:
		with_gen.to_csv('08_08_17_09_networks/genus/bins/'+gen_b_ndat[net], sep = '\t', index = False)

for net in range(len(spec_p_edges)):
	edges = pd.read_csv('08_08_17_09_networks/species/pears/'+spec_p_edges[net], sep = '\t')
	node_dat = pd.read_csv('08_08_17_09_networks/species/pears/'+spec_p_ndat[net], sep = '\t')
	with_gen = make_meta_from_file(edges, gender, spec_lv, existing_table = node_dat)
	if len(with_gen)> 0:
		with_gen.to_csv('08_08_17_09_networks/species/pears/'+spec_p_ndat[net], sep = '\t', index = False)

for net in range(len(spec_b_edges)):
	edges = pd.read_csv('08_08_17_09_networks/species/bins/'+spec_b_edges[net], sep = '\t')
	node_dat = pd.read_csv('08_08_17_09_networks/species/bins/'+spec_b_ndat[net], sep = '\t')
	with_gen = make_meta_from_file(edges, gender, spec_lv, existing_table = node_dat)
	if len(with_gen)> 0:
		with_gen.to_csv('08_08_17_09_networks/species/bins/'+spec_b_ndat[net], sep = '\t', index = False)

