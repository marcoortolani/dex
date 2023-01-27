# https://www.statology.org/seaborn-stacked-bar-plot/

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

cmpfile = '../data/Gene_biotype_summary.tsv'
geneCMP = pd.read_csv(cmpfile, index_col=0, sep='\t', comment='#')
geneCMP = geneCMP.transpose()

print(geneCMP.head())

# set seaborn plotting aesthetics
sns.set(style='white')

# Set your custom color palette
customPalette = sns.set_palette(sns.color_palette(n_colors=geneCMP.shape[1]))

#create stacked bar chart
geneCMP[0:2].plot(kind='bar', stacked=True, palette=customPalette)

#add overall title
plt.title('Customers by Time & Day of Week', fontsize=16)

#add axis titles
plt.xlabel('Day of Week')
plt.ylabel('Number of Customers')

#rotate x-axis labels
plt.xticks(rotation=45)

plt.show()