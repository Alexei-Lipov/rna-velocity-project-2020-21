# scVelo20.py but plots heatmaps for top 50 genes for each cluster.

# scVelo20.py but with higher louvain resolution (trying to get more than 9 clusters), and some lines commented out (plots and writing h5ad file)

# scVelo11.py but with the new larger dataset 18/02/21 (169k cells)

# VERSION WITH CLUSTER LABELLING & UNTRUNCATED DATASET & WITH SCANPY PREPROCESSING RATHER THAN SCTRANSFORM # WITH GENE ANALYSIS AND SOME LOUVAIN CLUSTERING AND GENE ANALYSIS ON THESE CLUSTERS, ALSO WRITES A H5AD FILE AFTER CALCULATING VELOCITY / VELOCITY GRAPH SO FOR FUTURE SCRIPTS CAN JUST LOAD IT IN - COULD EVEN DO IT VIA SH PROBABLY. ALSO PROPORTIONS PLOT.

#pip install -U scvelo --user
#pip install -U scanpy --user
#pip install git+https://github.com/theislab/scvelo --user
#pip install -U hnswlib --user
#pip install -U python-igraph --upgrade --quiet --user
#pip install -U python-igraph louvain --user



import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os

print("current working directory:")
print(os.getcwd())


scv.set_figure_params()

#adata = scv.read(filename = "/data/phar-ta-heart/lina3770/data/Data/unspliced_unprocessed_untruncated.csv", cache=True)
adata = scv.read(filename = "/data/phar-ta-heart/lina3770/data/Data/unspliced_unprocessed_untruncated_150221.csv", cache=True)

adata = adata.T

adata.layers["unspliced"] = adata.X


#adata.write(filename = "untruncated_adata_test.h5ad")

#adata = scv.read(filename = "/data/phar-ta-heart/lina3770/untruncated_adata_test.h5ad", cache=True)

#ldata = scv.read(filename = "/data/phar-ta-heart/lina3770/data/Data/spliced_unprocessed_untruncated.csv", cache=True)
ldata = scv.read(filename = "/data/phar-ta-heart/lina3770/data/Data/spliced_unprocessed_untruncated_150221.csv", cache=True)
ldata = ldata.T
ldata.layers["spliced"] = ldata.X

print("merge")
adata = scv.utils.merge(adata, ldata)

#adata = scv.datasets.pancreas()

#adata = adata[1:6000,1:500]

#3000 / n neighbours 400 works

#6000 / n neighbours 800 works
print("adata")
print(adata) 
print("adata.layers")
print(adata.layers)
print("adata.obs")
print(adata.obs)


#adata.obs.index = adata.obs["sample_batch"].astype(str) + adata.obs.index.astype(str)

#print(adata.obs)

#df['col'] = 'str' + df['col'].astype(str)



#AMs = pd.read_csv("/data/phar-ta-heart/lina3770/AM_names.csv")
#VMs = pd.read_csv("/data/phar-ta-heart/lina3770/VM_names.csv")

#AMs["Cluster"] = ["Atrial Myocyte"] * len(AMs)
#VMs["Cluster"] = ["Ventricular Myocyte"] * len(VMs)

#clusters = pd.concat([AMs, VMs])

#del clusters['Unnamed: 0']

#clusters = clusters.rename(columns={'x':"Cell"})

#print(clusters)

#print(len(adata.obs.index))
#print(len(clusters["Cell"]))

#intersection = list(set(adata.obs.index) & set(clusters["Cell"]))

#print(intersection)
#print(len(intersection))

#clusters = clusters[clusters["Cell"].isin(intersection)]

#print(len(clusters))
#print(clusters)

#print(clusters.index)

#clusters = clusters.set_index("Cell")

#print(clusters)

#print(clusters.index)

#clusters = clusters.reindex(adata.obs.index)

#print(clusters)

#adata.obs["clusters"] = list(clusters["Cluster"])

#print(adata.obs)

print("Debug")
print(adata.layers["spliced"][1:80,1:80])
print(adata.layers["spliced"].shape)
print(adata.layers["unspliced"][1:80,1:80])
print(adata.layers["unspliced"].shape)

print("amax")
print(np.amax(adata.layers["unspliced"],0))

scv.utils.show_proportions(adata)

scv.pp.filter_and_normalize(adata)

scv.utils.show_proportions(adata)

#scv.pl.proportions(adata, show = False)

print("proportions")
#plt.savefig("proportions_piechart.pdf")

#adata = adata[:,]

scv.pp.remove_duplicate_cells(adata)

#sc.pp.highly_variable_genes(adata)


#scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
# try normalizing in scanpy



scv.pp.moments(adata, use_highly_variable=False)

print(adata)

print(adata.layers["Ms"])
print(adata.layers["Mu"])
print(adata.X.shape)




scv.tl.velocity(adata)


pca = PCA(n_components=2)
reduced = pca.fit_transform(adata.layers["unspliced"])
t = reduced.transpose()
print("debug 4")
plt.scatter(t[0], t[1], s = 0.001)
plt.savefig("pca_rna_vel_untrunc_unspliced_169k.pdf")

plt.close() 

pca = PCA(n_components=2)
reduced = pca.fit_transform(adata.layers["spliced"])
t = reduced.transpose()
print("debug 4")
plt.scatter(t[0], t[1], s = 0.001)
plt.savefig("pca_rna_vel_untrunc_spliced_169k.pdf")



#scv.pl.scatter(adata.layers["velocity"], save = "pl_scatter.pdf")

#sc.pl.pca(adata, save = "pca_truncated_dataset.pdf")



print(adata)
print(adata.layers["velocity"])
print(adata.layers["velocity"].shape)


#pca = PCA(n_components=2)
#reduced = pca.fit_transform(adata.layers["velocity"])
#t = reduced.transpose()
print("debug 4")
#plt.scatter(t[0], t[1], s = 0.001)
#plt.savefig("pca_rna_vel_trunc.pdf")
print("debug 3")
#sc.tl.pca(adata.layers["velocity"])

print(adata)
print("debug 2")
#sc.pl.pca(adata.layers["velocity"], save = "pca_rna_vel_truncated_dataset.pdf")

print("debug 1")

#scv.pl.velocity(adata, save = "truncated_plt_RNA_vel.pdf")
scv.pp.neighbors(adata, n_neighbors = 800, method = 'hnsw')
#, n_neighbors = 100, knn = True,


print(adata)

print(adata.layers["Ms"])
print(adata.layers["Mu"])
print(adata.X.shape)


scv.tl.velocity_graph(adata)

#scv.pl.velocity_graph(adata, save = "velocity_graph_scanpy_processed.pdf")

scv.tl.umap(adata)

print(adata)

print(adata.obs)

#print("writing untruncated_adata_scanpy_and_scvelo_processed_deterministic_169k.h5ad")
#adata.write(filename = "untruncated_adata_scanpy_and_scvelo_processed_deterministic_169k.h5ad")




scv.tl.louvain(adata, resolution=4.0)            #resolution = 2.0 works and gives 16 clusters

print(adata)

print(adata.obs)

print("writing louvain clustering")
pd.set_option("display.max_rows", None, "display.max_columns", None)

f = open("169k_Louvain_Clustering_higher_res.txt", "w")
f.write(str(adata.obs["louvain"]))
f.close()


sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save = "untruncated_RNA_vel_scanpy_proc_gene_rankings_louvain_clusters_169k_res4_scvelo21.pdf")


scv.tl.rank_velocity_genes(adata, groupby='louvain', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()



# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='louvain')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.05,
            min_edge_width=2, node_size_scale=0.5, save = "untruncated_RNA_vel_scanpy_proc_paga_169k_res4_scvelo21.pdf",
            legend_fontsize = 8, legend_fontweight = 'normal', legend_loc = "on data")

scv.pl.velocity_embedding(adata, basis = "umap", arrow_length=0.7, color = ["louvain"], arrow_size=0.5, dpi=1000, fontsize = 1, title = "", save = "untruncated_RNA_vel_clustered_scanpy_processed_louvain_169k_res4_scvelo21.png", legend_fontsize = 8, legend_fontweight = 'normal', legend_loc = "on data")
#scv.pl.velocity_embedding(adata, basis = "umap", arrow_length=0.7, arrow_size=0.5, legend_loc = "right margin", dpi=1000, fontsize = 5, title = "", save = "untruncated_RNA_vel_clustered_scanpy_processed.pdf")
#scv.pl.velocity_embedding_grid(adata, basis='umap', save = "untruncated_RNA_vel_grid_clustered_scanpy_processed.pdf")
scv.pl.velocity_embedding_stream(adata, basis='umap', color = ["louvain"], min_mass = 0, fontsize = 1, linewidth = 0.3, size = 1, dpi=1000, save = "untruncated_RNA_vel_stream_clustered_scanpy_processed_by_louvain_169k_res4_scvelo21.pdf", legend_fontsize = 8, legend_fontweight = 'normal', legend_loc = "on data")

num_clusters = len(adata.obs["louvain"].cat.categories)
print(num_clusters)

for i in range(0,num_clusters):
  # SLICING TO GET ONLY CLUSTER i CELLS FROM ADATA BASED ON CLUSTER CATEGORY "i"
  Heatmap_adata = adata[adata.obs["louvain"] == adata.obs["louvain"].cat.categories[i]]
  sc.pl.heatmap(Heatmap_adata, adata.uns['rank_velocity_genes']['names'][str(i)][0:50], groupby='louvain', cmap='viridis', dendrogram=False, save = "untruncated_RNA_vel_scanpy_proc_169k_res4_heatmap_cluster" + str(i) + ".pdf")
  
  
print("writing untruncated_adata_scanpy_and_scvelo_processed_deterministic_169k_louvain_clustered_res4.h5ad")
adata.write(filename = "untruncated_adata_scanpy_and_scvelo_processed_deterministic_169k_louvain_clustered_res4.h5ad")