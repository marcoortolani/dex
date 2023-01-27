import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np

# Reading in data from "new" csv file (the "Feature ID" column will be used as id)
dexfile = '../data/Merged_differential_expression_results_CMM.csv'
dex = pd.read_csv(dexfile, index_col=0, comment='#')

# https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html

plt.scatter(x=dex['BL21_3hr-control_3hr_logFC'],y=dex['BL21_3hr-control_3hr_adj.P.Val'].apply(lambda x:-np.log10(x)),s=1)



plt.title("BL21 vs control 3h")
plt.xlabel("log2FC")
plt.ylabel("adjPVal")
plt.axvline(-2,color="grey",linestyle="--")
plt.axvline(2,color="grey",linestyle="--")
plt.axhline(2,color="grey",linestyle="--")

#plt.show()

plt.savefig('../plots/volcano.png')
