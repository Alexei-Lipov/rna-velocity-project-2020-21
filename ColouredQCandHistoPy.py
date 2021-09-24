import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white', fontsize=10, figsize=[8, 8])

print("loading")
adata_1 = sc.read_csv('/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E8_5_1/E8_5_1_Spliced.csv') 
adata_1 = adata_1.T
print("loading")
adata_2 = sc.read_csv('/data/phar-ta-heart/lina3770/data/Data/Spliced_Unspliced/E8_5_2/E8_5_2_Spliced.csv')
adata_2 = adata_2.T

#adata_1.var_names_make_unique()
#t = adata_1.var_names.str.endswith('-1')
#for i in range(1,16):
#    t = np.add(t, adata_1.var_names.str.endswith('-'+str(i+1)))
#t = ~t
#adata_1 = adata_1[:, t]

for j in ["adata_1", "adata_2"]:
  j.var_names_make_unique()
  t = j.var_names.str.endswith('-1')
  for i in range(1,16):
      t = np.add(t, j.var_names.str.endswith('-'+str(i+1)))
  t = ~t
  j = j[:, t]

adata = adata_1.concatenate(adata_2)

adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Number of genes expressed by each cell versus total counts for each cell. The colour of the cell corresponds to the fraction of its counts that are mitochondrial, with lighter colours indicating a higher percentage.
sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_mt', title = '')

# Zooming into the bottom-left corner:
sc.pl.scatter(adata[adata.obs.n_genes_by_counts < 600, :], 'total_counts', 'n_genes_by_counts', color='pct_counts_mt', title = '')
