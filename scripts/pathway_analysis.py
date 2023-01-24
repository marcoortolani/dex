###########
# Analysis of pathways for the experiments on DESCR
# and production of plots:
# - pair plots 
# - volcano plot
#
# Conda env: celeste
# MO - Jan 2023
###########

# Import necessary Python libraries

# Data analysis
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

# Reading in data from "new" csv file (the "Feature ID" column will be used as id)
dexfile = '../data/Merged_differential_expression_results_CMM.csv'
dex = pd.read_csv(dexfile, index_col=0, comment='#')

KeyPathway_IL_17 = ('CCL2', 'CCL20', 'CSF2', 'CSF3', 'CXCL1', 
                    'CXCL10', 'CXCL2', 'CXCL3', 'CXCL8', 'IL1B', 
                    'IL6', 'MMP3', 'PTGS2', 'TNF', 'TNFAIP3')

KeyPathway_TNF = ('CCL2', 'CCL20', 'CSF2', 'CX3CL1', 'CXCL1', 'CXCL10', 
                'CXCL2', 'CXCL3', 'IFNB1', 'IL1B', 'IL6', 'IRF1', 'LIF', 
                'MMP3', 'NOD2', 'PTGS2', 'SELE', 'TNF', 'TNFAIP3', 'VCAM1')

KeyPathway_TollLike = ('CCL3', 'CCL4', 'CXCL10', 'CXCL11', 
                        'CXCL8', 'IFNB1', 'IL1B', 'IL6', 'TNF')


colnames = pd.DataFrame(columns=['HRPII_BL21', 'HRPII_CC', 'LPS'])
colnames.loc['logFC_3h'] = ['BL21_3hr-control_3hr_logFC', 'CC_3hr-control_3hr_logFC', 'LPS_3hr-control_3hr_logFC']
colnames.loc['logFC_12h'] = ['BL21_12hr-control_12hr_logFC', 'CC_12hr-control_12hr_logFC', 'LPS_12hr-control_12hr_logFC']
colnames.loc['adjPval_3h'] = ['BL21_3hr-control_3hr_adj.P.Val', 'CC_3hr-control_3hr_adj.P.Val', 'LPS_3hr-control_3hr_adj.P.Val'] 
colnames.loc['adjPval_12h'] =['BL21_12hr-control_12hr_adj.P.Val', 'CC_12hr-control_12hr_adj.P.Val', 'LPS_12hr-control_12hr_adj.P.Val']
colnames.loc['sample1_12h'] = ['sample.24', 'sample.16', 'sample.32']
colnames.loc['sample2_12h'] = ['sample.23', 'sample.15', 'sample.31']
colnames.loc['sample3_12h'] = ['sample.22', 'sample.14', 'sample.30']
colnames.loc['sample4_12h'] = ['sample.21', 'sample.13', 'sample.29']
colnames.loc['sample1_3h'] = ['sample.20', 'sample.12', 'sample.28']
colnames.loc['sample2_3h'] = ['sample.19', 'sample.11', 'sample.27']
colnames.loc['sample3_3h'] = ['sample.18', 'sample.10', 'sample.26']
colnames.loc['sample4_3h'] = ['sample.17', 'sample.9', 'sample.25']

for experiment in ['logFC_3h', 'logFC_12h']:
    selectedcols = colnames.loc[experiment]
    pathdata = dex[selectedcols]
    pathdata['external_gene_name'] = dex['external_gene_name']

    # pathdata = pathdata[pathdata['external_gene_name'].isin(KeyPathway_IL_17)]

    g = sns.pairplot(pathdata)
    plt.savefig('../plots/'+ experiment+'.png')



exit()
mde_3h_BL21_CC_LPS = pd.read_csv('../data/MDE_Log2FC_3h_BL21_CC_LPS.csv', index_col=0, comment='#')

# col_name = names.loc['3h_BL21_CC_LPS']['HRPII_BL21']

# for x in KeyPathway_IL_17_12h:
#     print(mde_3h_BL21_CC_LPS.loc[x][col_name])

# plt.figure(figsize=(9,6)) # Set plot dimensions

# sns.regplot(x=mde_3h_BL21_CC_LPS.index, y="BL21_3hr-control_3hr_logFC", data=mde_3h_BL21_CC_LPS)

#pathdata = mde_3h_BL21_CC_LPS[mde_3h_BL21_CC_LPS.index.isin(KeyPathway_IL_17_12h)]
pathdata = mde_3h_BL21_CC_LPS

pathdata = pathdata.reset_index(level=0)
melted_df = pd.melt(pathdata, id_vars="external_gene_name", var_name="experiment")

sns.catplot(data=melted_df, x="external_gene_name", y="value", hue="experiment", kind="point")

# https://seaborn.pydata.org/tutorial/categorical.html
# https://stats.stackexchange.com/questions/108007/correlations-with-unordered-categorical-variables
# https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/iris_plots/
# https://pythonbasics.org/seaborn-pairplot/

g = sns.pairplot(pathdata)
plt.savefig('save_as_a_png.png')


logFCcols = ['BL21_3hr-control_3hr_logFC', 'CC_3hr-control_3hr_logFC', 'LPS_3hr-control_3hr_logFC']

corr_matrix = pathdata[logFCcols].corr()

print(corr_matrix)

sns.heatmap(corr_matrix, annot=True)

plt.show()