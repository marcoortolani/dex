#!/usr/bin/env python
# coding: utf-8

# # Differential Expression Analysis
#
# This notebook is a first attempt to automate differential gene expression analysis with Python.
#
# Gene expression data are read from the .tsv file and imported into Python data structures so that we can select "interesting" genes to submit to Enrichr for enrichment analysis.
#
# We'll follow these steps:
# 1. Filter the most significant genes (p<0.01)
# 2. Filter Log2FC (Less/equal to -1 or greater/equal to 1)
# 3. Clean up data and maintain only Entrez and Ensemble ID, Log2 FC, p.value cols
# 4. Rank the genes (largest to  smallest using Log2FC)
# 5. Paste results on Enrichr

# In[ ]:


# Import necessary Python libraries

# Data analysis
import numpy as np
import pandas as pd

# Bioinformatics
import gseapy as gp

# Graphics and plotting
import matplotlib.pyplot as plt

#import seaborn as sns


# ## Data preparation

# In[ ]:


# Reading in data from tsv file (the "Feature ID" column will be used as id)

filename = '../data/Merged_differential_expression_results.tsv'

#dex = pd.read_csv(filename, index_col=0, sep='\t')
dex = pd.read_csv(filename, dtype={'entrezgene': str}, sep='\t')


# In[ ]:


# Separate the columns into strings (all gene names, descriptions, ...) and numeric (gene expressions, fold changes, and quality metrics)
all_cols = list(dex.columns)
descr_cols = ['Feature_ID', 'entrezgene', 'external_gene_name', 'gene_biotype', 'external_gene_source', 'description']
numeric_cols = [x for x in all_cols if x not in set(descr_cols)]

# The select samples, fold metrics
sample_cols = [x for x in numeric_cols if 'sample' in x]
linfold_cols = [x for x in numeric_cols if 'linearFC' in x]
logfold_cols = [x for x in numeric_cols if 'logFC' in x]

control_3h_cols = ['sample.'+str(x) for x in range(1,5)]
control_12h_cols = ['sample.'+str(x) for x in range(5,9)]
CC_3h_cols = ['sample.'+str(x) for x in range(9,13)]
CC_12h_cols = ['sample.'+str(x) for x in range(13,17)]
BL21_3h_cols = ['sample.'+str(x) for x in range(17,21)]
BL21_3h_bcols = ['sample.'+str(x) for x in range(21,25)]
LPS_3h_cols = ['sample.'+str(x) for x in range(25,29)]
LPS_21h_cols = ['sample.'+str(x) for x in range(29,33)]

# Associate fold column name with p-value column name
adj_p_value_columns = {
    'BL21_12hr-BL21_3hr_logFC': 'BL21_12hr-BL21_3hr_adj.P.Val',
    'BL21_12hr-control_12hr_logFC': 'BL21_12hr-BL21_3hr_adj.P.Val',
    'BL21_3hr-control_3hr_logFC': 'BL21_3hr-control_3hr_adj.P.Val',
    'CC_12hr-CC_3hr_logFC': 'CC_12hr-CC_3hr_adj.P.Val',
    'CC_12hr-control_12hr_logFC':'CC_12hr-control_12hr_adj.P.Val',
    'CC_3hr-control_3hr_logFC':'CC_3hr-control_3hr_adj.P.Val',
    'control_12hr-control_3hr_logFC':'control_12hr-control_3hr_adj.P.Val',
    'LPS_12hr-control_12hr_logFC':'LPS_12hr-control_12hr_adj.P.Val',
    'LPS_12hr-LPS_3hr_logFC':'LPS_12hr-LPS_3hr_adj.P.Val',
    'LPS_3hr-control_3hr_logFC':'LPS_3hr-control_3hr_adj.P.Val'
}


# ### Producing ranking files for all experiments data

# In[ ]:


# Thresholds for selecting the relevant rows

# CHANGE VALUES HERE AS NEEDED (then re-run: Cell -> Run all)
pValThreshold = 0.05
logFCThreshold = 2

# Uncomment the line corresponding to the gene name you'd like in the output (and comment the others)
#geneCol = 'external_gene_name'
#geneCol = 'Feature_ID'
#geneCol = 'entrezgene'
geneCol = 'external_gene_name'


# In[ ]:


print("Selecting only significant experiments with adj p-val < " + str(pValThreshold) + " - |logFC| <= " +str(logFCThreshold))
for data_column_name in logfold_cols:
    pval_column_name = adj_p_value_columns[data_column_name]


    # Select only significant p-values (rows for which p-value is less than a threshold, e.g.: p<0.01)
    data1 = dex[dex[pval_column_name]<pValThreshold]

    count1 = data1.shape[0]

    # Select only LogFC values which are <=-logFCThreshold or >=logFCThreshold (e.g.: |logFC| <= 1 )
    data2 = data1[(data1[data_column_name]<=-logFCThreshold) | (data1[data_column_name]>=logFCThreshold)]

    count2 = data2.shape[0]

    # print(f'{data_column_name:<30} : {str(count1):>6} (after p-val) - {str(count2):>6} (after LogFC)')

    # get rid of all columns we are not interested in right now
    cols_to_save = [geneCol, data_column_name, pval_column_name]
    col_labels = ['GeneID', 'Log2FC', 'AdjPval']

    data1 = data1[cols_to_save]
    data2 = data2[cols_to_save]

    # Now sort rows by ascending expression value
    data1 = data1.sort_values(
        by=data_column_name,
        ascending=False
        )

    #data1.to_csv('../results/'+data_column_name+'-v1.rnk', index=False, sep='\t', header=False)
    #data2.to_csv('../results/'+data_column_name+'-v2.rnk', index=False, sep='\t', header=False)

    data2.to_csv('../results/'+data_column_name+'.tsv', index =False, sep='\t', header=col_labels)

print()
#print("Saved all files: v1 is with pVal threshold only; v2 with pVal and LogFC thresholds")
print("Saved all files: check the results directory")


# You can find all results in the `data` folder (accessible via `File -> Open` then `../results`)

# In[ ]:
