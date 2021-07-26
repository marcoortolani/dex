#%matplotlib inline
#%config InlineBackend.figure_format='retina' # mac
#%load_ext autoreload
#%autoreload 2

import numpy as np
import pandas as pd
import seaborn as sns

import gseapy as gp
import matplotlib.pyplot as plt

filename = '../data/Merged_differential_expression_results.tsv'

#dex = pd.read_csv(filename, index_col=0, dtype={'entrezgene': 'str'}, sep='\t')
dex = pd.read_csv(filename, index_col=0, dtype={'entrezgene': 'str'}, sep='\t')

# Separate the columns into strings (all gene names, descriptions, ...) and numeric (gene expressions, fold changes, and quality metrics)
all_cols = list(dex.columns)
descr_cols = ['Feature_ID', 'entrezgene', 'external_gene_name', 'gene_biotype', 'external_gene_source', 'description']
numeric_cols = [x for x in all_cols if x not in set(descr_cols)]

# The select samples, fold metrics
sample_cols = [x for x in numeric_cols if 'sample' in x]
linfold_cols = [x for x in numeric_cols if 'linearFC' in x]
logfold_cols = [x for x in numeric_cols if 'logFC' in x]

control_3h_cols = ['sample.' + str(x) for x in range(1, 5)]
control_12h_cols = ['sample.' + str(x) for x in range(5, 9)]
CC_3h_cols = ['sample.' + str(x) for x in range(9, 13)]
CC_12h_cols = ['sample.' + str(x) for x in range(13, 17)]
BL21_3h_cols = ['sample.' + str(x) for x in range(17, 21)]
BL21_3h_bcols = ['sample.' + str(x) for x in range(21, 25)]
LPS_3h_cols = ['sample.' + str(x) for x in range(25, 29)]
LPS_21h_cols = ['sample.' + str(x) for x in range(29, 33)]

# Associate fold column name with p-value column name
adj_p_value_columns = {
    'BL21_12hr-BL21_3hr_logFC': 'BL21_12hr-BL21_3hr_adj.P.Val',
    'BL21_12hr-control_12hr_logFC': 'BL21_12hr-BL21_3hr_adj.P.Val',
    'BL21_3hr-control_3hr_logFC': 'BL21_3hr-control_3hr_adj.P.Val',
    'CC_12hr-CC_3hr_logFC': 'CC_12hr-CC_3hr_adj.P.Val',
    'CC_12hr-control_12hr_logFC': 'CC_12hr-control_12hr_adj.P.Val',
    'CC_3hr-control_3hr_logFC': 'CC_3hr-control_3hr_adj.P.Val',
    'control_12hr-control_3hr_logFC': 'control_12hr-control_3hr_adj.P.Val',
    'LPS_12hr-control_12hr_logFC': 'LPS_12hr-control_12hr_adj.P.Val',
    'LPS_12hr-LPS_3hr_logFC': 'LPS_12hr-LPS_3hr_adj.P.Val',
    'LPS_3hr-control_3hr_logFC': 'LPS_3hr-control_3hr_adj.P.Val'
}


data_column_name = 'control_12hr-control_3hr_logFC'
pval_column_name = adj_p_value_columns[data_column_name]


# Select only significant p-values
data = dex[dex[data_column_name]<0.01]

data = data[(data[data_column_name]<=-1) | (data[data_column_name]>=1)]

data = data[['entrezgene', 'external_gene_name', data_column_name]]

data = data.sort_values(
    by=data_column_name,
    ascending=True
    )

genes = data.head(10)[['entrezgene']]

genes.dropna(subset = ["entrezgene"], inplace=True)

#glist = genes.squeeze().str.strip().tolist()

names = gp.get_library_name() # default: Human
names[:10]

gene_list = pd.read_csv("../data/gseapy_gene_list.txt",header=None, sep="\t")
glist = gene_list.squeeze().str.strip().tolist()

enr = gp.enrichr(gene_list=glist,
                 gene_sets=['KEGG_2016','KEGG_2013'],
                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                 description='test_name',
                 outdir='test/enrichr_kegg',
                 #no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                 )

enr.results.head(5)