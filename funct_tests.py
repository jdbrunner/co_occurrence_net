from pylab import *
import pandas as pd
import sys
import itertools as iter
from scipy.special import binom as choose
from scipy.stats import binom
from numpy.random import binomial as bino
import re
import time
import numpy.ma as mask
from co_occ_funs import *


csv_name = 'merged_assignment.txt'

######Import the abundance matrix
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 
abundance_array_full = pd.read_csv(csv_name, sep = ' ')
#

#Combine samples of the same kind (that is, the same part of the body) so that we can 
#color according to where abundance of a genome is highest.
diff_samps_types = unique([name[:-10] for name in array(abundance_array_full.keys())[2:]])
data_samps_types = transpose([sum([abundance_array_full[samp] for samp in abundance_array_full.keys()[2:]
							 if smp_type in samp],axis = 0) for smp_type in diff_samps_types])

samp_type_abund = pd.DataFrame(data_samps_types, columns = diff_samps_types)

samp_type_norm = True
if samp_type_norm:
	for col in samp_type_abund.columns:
		tot = sum(samp_type_abund[col])
		samp_type_abund[col] = samp_type_abund[col]/tot


samp_type_abund['TAXA'] = abundance_array_full['TAXA']

s_type1 = color_picker(samp_type_abund.iloc[15][:-1])
print(s_type1)
