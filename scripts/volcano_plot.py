###########
# Producing volcano plots for the analysis of pathways.
# Either specific genes can be selected, or all genes
#
# Conda env: celeste
# MO - Jan 2023
###########

# Import necessary Python libraries
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

# Our thresholds of interest for p-value and fold count
pvalThres = 0.05
FCThres = 2

# Unselect for HRPII_BL21 (3h)
# logFCcol = colnames['HRPII_BL21']['logFC_3h']
# adjPvalcol = colnames['HRPII_BL21']['adjPval_3h']
# title = 'HRPII_BL21 (3h)'

# Unselect for HRPII_CC (3h)
# logFCcol = colnames['HRPII_CC']['logFC_3h']
# adjPvalcol = colnames['HRPII_CC']['adjPval_3h']
# title = 'HRPII_CC (3h)'

# Unselect for HRPII_BL21 (12h)
# logFCcol = colnames['HRPII_BL21']['logFC_12h']
# adjPvalcol = colnames['HRPII_BL21']['adjPval_12h']
# title = 'HRPII_BL21 (12h)'

# Unselect for HRPII_CC (12h)
# logFCcol = colnames['HRPII_CC']['logFC_12h']
# adjPvalcol = colnames['HRPII_CC']['adjPval_12h']
# title = 'HRPII_CC (12h)'

# Unselect for LPS (3h)
# logFCcol = colnames['LPS']['logFC_3h']
# adjPvalcol = colnames['LPS']['adjPval_3h']
# title = 'LPS (3h)'

# Unselect for LPS (12h)
logFCcol = colnames['LPS']['logFC_12h']
adjPvalcol = colnames['LPS']['adjPval_12h']
title = 'LPS (12h)'

# scatter plot for the selected experiment - non-significant genes
plt.scatter(x=dex[logFCcol],
            y=dex[adjPvalcol].apply(lambda x:-np.log10(x)),
            s=1,
            label="Not significant",
            color='grey')

# highlight down- or up- regulated genes
down = dex[(dex[logFCcol]<=-FCThres) & (dex[adjPvalcol]<=pvalThres)]
up = dex[(dex[logFCcol]>=FCThres) & (dex[adjPvalcol]<=pvalThres)]

# scatter plot for the selected experiment - down-regulated genes
plt.scatter(x=down[logFCcol],
            y=down[adjPvalcol].apply(lambda x:-np.log10(x)),
            s=3,
            label="Down-regulated",
            color="blue")

# scatter plot for the selected experiment - up-regulated genes
plt.scatter(x=up[logFCcol],
            y=up[adjPvalcol].apply(lambda x:-np.log10(x)),
            s=3,
            label="Up-regulated",
            color="red")

plt.title(title)
plt.xlabel("log$_2$FC")
plt.ylabel("-log$_{10}$adjPVal")

# adding 2 vertical lines to identify representative FC
plt.axvline(-FCThres,color="grey",linestyle="--")
plt.axvline(FCThres,color="grey",linestyle="--")

# adding horizontal line (and label) to mark significant p-value
ax = plt.gca()
yt = ax.get_yticks()
yt=np.append(yt,1.5)

ytl=yt.tolist()
ytl[-1]="-log$_{10}$(" + str(pvalThres) + ")"
ax.set_yticks(yt)
ax.set_yticklabels(ytl)

# plt.axhline(2,color="grey",linestyle="--")
plt.axhline(-np.log10(pvalThres),color="grey",linestyle="--")

plt.legend()

#plt.savefig('../plots/volcano_' + title+ '.png', bbox_inches="tight")
plt.savefig('../plots/volcano_' + title+ '.svg', bbox_inches="tight")
# plt.show()