# Import necessary Python libraries

# Data analysis
import numpy as np
import pandas as pd

# Bioinformatics
import gseapy as gp

# Graphics and plotting
import matplotlib.pyplot as plt
import seaborn as sns

# Reading in data from tsv file (the "Feature ID" column will be used as id)
filename = '../data/Merged_differential_expression_results.tsv'
dex = pd.read_csv(filename, index_col=0, sep='\t')

# Name all sample columns separately
control_3h_cols = ['sample.'+str(x) for x in range(1,5)]
control_12h_cols = ['sample.'+str(x) for x in range(5,9)]
CC_3h_cols = ['sample.'+str(x) for x in range(9,13)]
CC_12h_cols = ['sample.'+str(x) for x in range(13,17)]
BL21_3h_cols = ['sample.'+str(x) for x in range(17,21)]
BL21_12h_cols = ['sample.'+str(x) for x in range(21,25)]
LPS_3h_cols = ['sample.'+str(x) for x in range(25,29)]
LPS_12h_cols = ['sample.'+str(x) for x in range(29,33)]

data= pd.DataFrame(index=dex.index)

sel = dex[control_3h_cols]
data['control_3h'] = sel.mean(axis=1)
data['stddev'] = sel.std(axis=1)
# sel = dex[control_12h_cols]
# data['control_12h'] = sel.mean(axis=1)
# sel = dex[CC_3h_cols]
# data['CC_3h'] = sel.mean(axis=1)
# sel = dex[CC_12h_cols]
# data['CC_12h'] = sel.mean(axis=1)
# sel = dex[BL21_3h_cols]
# data['BL21_3h'] = sel.mean(axis=1)
# sel = dex[BL21_12h_cols]
# data['BL21_12h'] = sel.mean(axis=1)
# sel = dex[LPS_3h_cols]
# data['LPS_3h'] = sel.mean(axis=1)
# sel = dex[LPS_12h_cols]
# data['LPS_12h'] = sel.mean(axis=1)

data1 = dex[control_3h_cols]
data2 = dex[control_12h_cols]

sns.boxplot(data=data1)
plt.show()

xxx = data1.iloc[0]
print(type(xxx))

#sns.boxplot(data=xxx)
#plt.show()
