###########
# Plotting CPM per biotype as a stacked bar plot
#
# Conda env: celeste
# MO - Jan 2023
###########

# import libraries
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# Read in data from file 
cmpfile = '../data/Gene_biotype_summary.tsv'
geneCMP = pd.read_csv(cmpfile, index_col=0, sep='\t', comment='#')
geneCMP = geneCMP.transpose()

# reverse row order (just for clean plotting)
geneCMP = geneCMP.iloc[::-1]

# Choose an appropriate color palette
# I'll mix two palettes, to make it similar to our reference image
colors = [(0.97, 0.46, 0.42), (0, 0.74, 0.41)] 
top = LinearSegmentedColormap.from_list(
        "Custom", colors, N=19)
#colors = [(0, 0.74, 0.8), (0.99, 0.43, 0.52)] 
colors = [(0, 0.74, 0.498), (0.99, 0.43, 0.52)] 
bottom = LinearSegmentedColormap.from_list(
        "Custom", colors, N=18)

newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

# producing the plot
geneCMP.plot(kind = 'bar', stacked = True, width=0.8, colormap=newcmp)

# after plotting the data, format the labels
plt.rcParams['xtick.labelsize'] = 8
ticks_loc = plt.gca().get_yticks().tolist()
plt.gca().yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))

# using format string '{:.0f}' here but you can choose others
plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in ticks_loc])

# rename lables for x axis 
#plt.gca().set_xticklabels( [str(i) for i in reversed(range(1,33))] )
plt.gca().set_xticklabels( [str(i) for i in range(1,33)] )

#rotate x-axis labels
plt.xticks(rotation=45, ha='right', rotation_mode='anchor')

plt.xlabel("Sample")
plt.ylabel("Total gene CMP per biotype")

plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left", ncol=2,
            handlelength=1, handleheight=1,
            title="Gene biotype",
            fontsize="xx-small",
            frameon=False)

# plt.savefig("../plots/cpm_biotype.png", bbox_inches="tight")
plt.savefig("../plots/cpm_biotype.svg", bbox_inches="tight")

# show the graph
# plt.show()