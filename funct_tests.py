from pylab import *
import pandas as pd
import sys
from co_occ_funs import *
import os




# os.environ['JOBLIB_START_METHOD'] = 'forkserver'
# 
# import joblib

if __name__ == '__main__':

	lvl_abundance_array = pd.DataFrame(rand(100,50))
	lvl_abundance_array['LEVEL'] = array([10]*len(lvl_abundance_array))
	lvl_abundance_array['TAXA'] = rand(len(lvl_abundance_array))


	pears = build_network(lvl_abundance_array, 'pearson', thr = False, list_too = False)
	stats = mc_network_stats(lvl_abundance_array.values, pears.values, sims = 100)
	print(stats)





