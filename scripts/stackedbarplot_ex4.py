# https://www.python-graph-gallery.com/stacked-and-percent-stacked-barplot

# import libraries
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap

cmpfile = '../data/Gene_biotype_summary.tsv'
geneCMP = pd.read_csv(cmpfile, index_col=0, sep='\t', comment='#')
#geneCMP = pd.read_csv(cmpfile, index_col=None, comment='#', dtype={'gene_biotype': 'str'})
geneCMP = geneCMP.transpose()
# geneCMP = geneCMP.transpose().reset_index()
# geneCMP = geneCMP.rename({'index': 'Sample'}, axis='columns')

# bar1 = sns.barplot(x="Sample",  y="3prime_overlapping_ncRNA", data=geneCMP, color='darkblue')
# bar2 = sns.barplot(x="Sample", y="bidirectional_promoter_lncRNA", data=geneCMP, color='lightblue')

# # add legend
# top_bar = mpatches.Patch(color='darkblue', label='smoker = No')
# bottom_bar = mpatches.Patch(color='lightblue', label='smoker = Yes')
# plt.legend(handles=[top_bar, bottom_bar])

plt.rcParams['xtick.labelsize'] = 8
#plt.rcParams["axes.formatter.useoffset"] = False
# ax = plt.gca()
# ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

colors = sns.color_palette("viridis_r", n_colors=geneCMP.shape[1])
# colors = sns.color_palette("RdYlGn", n_colors=geneCMP.shape[1])
# colors = sns.color_palette("PiYG", n_colors=geneCMP.shape[1])
cmap1 = LinearSegmentedColormap.from_list("my_colormap", colors)

# geneCMP[['3prime_overlapping_ncRNA', 'bidirectional_promoter_lncRNA']].plot (kind = 'bar', stacked = True, color = ['darkblue', 'lightblue'])

# reverse row order (just for clean plotting)
geneCMP = geneCMP.iloc[::-1]

geneCMP.plot (kind = 'bar', stacked = True, width=0.8, colormap=cmap1)

# after plotting the data, format the labels
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

# plt.legend(loc='upper left')
#plt.legend(ncol=2)

# https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left", ncol=2,
            handlelength=1, handleheight=1,
            title="Gene biotype",
            fontsize="xx-small",
            frameon=False)

plt.savefig("../plots/cpm_biotype.png", bbox_inches="tight")

# show the graph
# plt.show()