# https://www.python-graph-gallery.com/stacked-and-percent-stacked-barplot

# import libraries
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
# from matplotlib.colors import DivergingNorm
from matplotlib.colors import TwoSlopeNorm
# checking default color
# sns.palplot(sns.color_palette(n_colors=37))

# diverging palette; two HUSL colors
#sns.palplot(sns.diverging_palette(250,10, 80, 50, center='dark', n=37))

# rdgn = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, sep=3, as_cmap=True)

# # normalize color
# vmin, vmax, vcenter = 0, 37, 20
# norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
# # create a normalized colorbar
# # cmap = 'RdBu'
# cbar = plt.cm.ScalarMappable(norm=norm, cmap=rdgn)

# sns.palplot(cbar)

top = mpl.colormaps['Oranges_r'].resampled(128)
#top = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, sep=3, as_cmap=True)
bottom = mpl.colormaps['Blues'].resampled(128)

colors = [(0.97, 0.46, 0.42), (0, 0.75, 0.73)] 
top = LinearSegmentedColormap.from_list(
        "Custom", colors, N=19)

colors = [(0, 0.74, 0.8), (0.99, 0.43, 0.52)] 
bottom = LinearSegmentedColormap.from_list(
        "Custom", colors, N=18)

newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

sns.palplot(newcmp(np.linspace(0,1,37)))
plt.show()