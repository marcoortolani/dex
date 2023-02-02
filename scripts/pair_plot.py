###########
# Producing pair plots for the analysis of pathways.
# Either specific genes can be selected, or all genes
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

colPlotnames = pd.DataFrame(columns=['HRPII_BL21', 'HRPII_CC', 'LPS'])
colPlotnames.loc['logFC_3h'] = ['HRPII_BL21 (3h)', 'HRPII_CC (3h)', 'LPS (3h)']
colPlotnames.loc['logFC_12h'] = ['HRPII_BL21 (12h)', 'HRPII_CC (12h)', 'LPS (12h)']
colPlotnames.loc['adjPval_3h'] = ['BL21_3hr-control_3hr_adj.P.Val', 'CC_3hr-control_3hr_adj.P.Val', 'LPS_3hr-control_3hr_adj.P.Val']
colPlotnames.loc['adjPval_12h'] =['BL21_12hr-control_12hr_adj.P.Val', 'CC_12hr-control_12hr_adj.P.Val', 'LPS_12hr-control_12hr_adj.P.Val']
colPlotnames.loc['sample1_12h'] = ['sample.24', 'sample.16', 'sample.32']
colPlotnames.loc['sample2_12h'] = ['sample.23', 'sample.15', 'sample.31']
colPlotnames.loc['sample3_12h'] = ['sample.22', 'sample.14', 'sample.30']
colPlotnames.loc['sample4_12h'] = ['sample.21', 'sample.13', 'sample.29']
colPlotnames.loc['sample1_3h'] = ['sample.20', 'sample.12', 'sample.28']
colPlotnames.loc['sample2_3h'] = ['sample.19', 'sample.11', 'sample.27']
colPlotnames.loc['sample3_3h'] = ['sample.18', 'sample.10', 'sample.26']
colPlotnames.loc['sample4_3h'] = ['sample.17', 'sample.9', 'sample.25']

for experiment in ['logFC_3h', 'logFC_12h']:
    selectedcols = colnames.loc[experiment].to_list()
    selectedcols.append('external_gene_name')
    pathdata = dex[selectedcols]
    
    # rename for plotting
    newColNames = colPlotnames.loc[experiment].to_list()
    newColNames.append('external_gene_name')
    pathdata.columns = newColNames
    # pathdata.set_axis(newColNames, axis=1, copy=False)

    # pathdata = pathdata[pathdata['external_gene_name'].isin(KeyPathway_IL_17)]

    pp = sns.pairplot(pathdata)

    # pp.axes[0,0].xlabel("log2FC")
#    plt.savefig('../plots/'+ experiment+'.png')
    plt.savefig('../plots/'+ experiment+'.svg')

