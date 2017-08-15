#########################################################################################
#																						#		
#				Co-Occurrence Network Color Key Builder				 					#
#																						#
#										Jim Brunner 									#
#																						#
#																						#
#																						#
#########################################################################################

#Quick script to display the colors used in co_occurrence.py

## modules needed
from pylab import *
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors


#####Variables
##User input, should be entered as an argument, should be same as used in co_occurrence.py
csv_name = sys.argv[1]

######Import the abundance matrix
#Makes a pandas DataFrame, use abundance_array.keys() to see column labels, 
abundance_array = pd.read_csv(csv_name, sep = ' ')[:1]

diff_samps_types = unique([name[:-10] for name in array(abundance_array.keys())[2:]])
numcols = len(diff_samps_types)
colors = linspace(0,255,numcols)

rgb_values = [cm.rainbow(col.astype(int)) for col in colors]
#hex_values = [matplotlib.colors.rgb2hex(rgb) for rgb in rgb_values]
the_cols = dict()
for i in diff_samps_types:
	the_cols[i] = rgb_values[argwhere(diff_samps_types == i)]

n = len(diff_samps_types)
ncols = 4
nrows = n // ncols + 1

fig, ax = plt.subplots(figsize=(8, 5))

# Get height and width
X, Y = fig.get_dpi() * fig.get_size_inches()
h = Y / (nrows + 1)
w = X / ncols

for i, name in enumerate(diff_samps_types):
    col = i % ncols
    row = i // ncols
    y = Y - (row * h) - h

    xi_line = w * (col + 0.05)
    xf_line = w * (col + 0.25)
    xi_text = w * (col + 0.3)

    ax.text(xi_text, y, name, fontsize=(h * 0.1),
            horizontalalignment='left',
            verticalalignment='center')

    ax.hlines(y + h * 0.1, xi_line, xf_line,
              color=the_cols[name], linewidth=(h * 0.6))

ax.set_xlim(0, X)
ax.set_ylim(0, Y)
ax.set_axis_off()

fig.subplots_adjust(left=0, right=1,
                    top=1, bottom=0,
                    hspace=0, wspace=0)
plt.show()