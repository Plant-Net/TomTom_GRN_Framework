import sys
import numpy as np
import gudhi as gd
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.metrics import pairwise_distances
from sklearn.manifold import MDS
from sklearn.preprocessing import LabelEncoder
from gudhi.clustering.tomato import Tomato
from sklearn.base import BaseEstimator, TransformerMixin
from scipy.sparse.csgraph import connected_components
from gudhi.cover_complex import MapperComplex
import pickle as pck
import networkx as nx
import seaborn as sns
from matplotlib.colors import Normalize
import math
from networkx.drawing.nx_agraph import graphviz_layout
plt.rcParams['svg.fonttype'] = 'none'
# Clustering method based on GRN
class GraphCC(BaseEstimator, TransformerMixin):
    def __init__(self, graph):
        self.graph = graph
    def fit(self, X, y=None):
        return self
    def fit_predict(self, X, y=None):
        n_components, labels = connected_components(csgraph=self.graph[:,X[:,0]][X[:,0],:], directed=False)
        self.n_clusters_ = n_components
        self.labels_ = labels
        return self.labels_

# Read GRN
X = pd.read_csv('./Curated_gene_regulatory_network.tsv', sep='\t', header=0, index_col=None)
X = np.array(X)
# print(len(X), len(np.unique(X[:,0])), len(np.unique(X[:,1])))

# Encoding names with integers
le = LabelEncoder()
n_tf, n_genes = len(np.unique(X[:,0])), len(np.unique(X[:,1]))
le.fit(np.array(list(np.unique(X[:,0])) + list(np.unique(X[:,1]))))

X[:,0] = le.transform(X[:,0])
X[:,1] = le.transform(X[:,1])

# Read expression + stats
F1 = pd.read_csv('./LogFC_genes_in_cGRN.tsv', sep='\t', header=0, index_col=None)
F2 = pd.read_csv('./Stats_genes_in_cGRN.tsv', sep='\t', header=0, index_col=None)
# print(np.array(F1), np.array(F2))
# print(F1.columns)

F1, F2 = np.array(F1), np.array(F2)
F1[:,0] = le.transform(F1[:,0])
F2[:,0] = le.transform(F2[:,0])
# print(F1,F2)

# Read activity
A = pd.read_csv('./Activity_TF_in_cGRN.tsv', sep='\t', header=0, index_col=None)
A = np.array(A)
A[:,0] = le.transform(A[:,0])

### TF-only GRN

# Restrict to TF
f1_tf,_,tf_indices = np.intersect1d(X[:,0], F1[:,0], return_indices=True)
F1_tf = F1[tf_indices,:]
f2_tf,_,tf_indices = np.intersect1d(X[:,0], F2[:,0], return_indices=True)
F2_tf = F2[tf_indices,:]
a_tf,_,tf_indices = np.intersect1d(X[:,0], A[:,0], return_indices=True)
A_tf = A[tf_indices,:]

le_tf = LabelEncoder()
le_tf.fit(X[:,0])

f1_perm = le_tf.transform(f1_tf)
inv_perm = np.array([np.argwhere(f1_perm == i)[0][0] for i in range(n_tf)]) 
F1_tf = F1_tf[inv_perm]
f2_perm = le_tf.transform(f2_tf)
inv_perm = np.array([np.argwhere(f2_perm == i)[0][0] for i in range(n_tf)]) 
F2_tf = F2_tf[inv_perm]
A_perm = le_tf.transform(a_tf)
inv_perm = np.array([np.argwhere(A_perm == i)[0][0] for i in range(n_tf)]) 
A_tf = A_tf[inv_perm]

common_tf = np.intersect1d(X[:,0], X[:,1])
targets_tf = [[] for _ in range(n_tf)]
graph_tf = np.zeros([n_tf, n_tf])
common_targets_tf = np.zeros([n_tf, n_tf])

# Construct networkx graph
for pair in X:
    tf0_id = le_tf.transform([pair[0]])[0]
    if pair[1] in common_tf:
        tf1_id = le_tf.transform([pair[1]])[0]
        graph_tf[tf0_id, tf1_id] = 1
        graph_tf[tf1_id, tf0_id] = 1
    targets_tf[tf0_id].append(pair[1])

# Compute common gene targets between TFs
for i in range(n_tf):
    for j in range(i+1,n_tf):
        common_targets_tf[i,j] = len(np.intersect1d(np.array(targets_tf[i]), np.array(targets_tf[j])))
        common_targets_tf[j,i] = len(np.intersect1d(np.array(targets_tf[i]), np.array(targets_tf[j])))

common_targets_tf_save = common_targets_tf

# plt.figure(figsize=(12, 10))  # Adjust the figure size for better visibility
# sns.heatmap(
#     common_targets_tf_save,
#     cmap="magma_r",  # Choose a colormap
#     annot=False,      # Disable annotations for large arrays
#     cbar=True,        # Show the colorbar
#     cbar_kws={'label': 'Value'}
# )


common_targets_tf = np.where(common_targets_tf > 0, np.ones(common_targets_tf.shape), np.zeros(common_targets_tf.shape))

#G = nx.Graph(graph_tf)
#plt.figure()
#nx.draw(G, node_size=1)
#plt.show()

# Initialize clustering based on GRNs
# Common targets-based clustering
grnCC_tf = GraphCC(common_targets_tf)
# GRN-based clustering
# grnCC_tf = GraphCC(graph_tf)

# Build filter out of graph Laplacian
n_targets_tf = np.array([len(trg) for trg in targets_tf])
L = nx.normalized_laplacian_matrix(nx.Graph(graph_tf))
eigenvalues, eigenvectors = np.linalg.eig(L.toarray())
eigenvectors = np.real(eigenvectors)
print(eigenvectors.shape)

#plt.figure()
#plt.hist(n_targets_tf, bins=300)
#plt.show()
#plt.figure()
#plt.hist(eigenvectors[:,0], bins=300)
#plt.show()


A_tf_filtered = np.where(np.abs(A_tf) < 2, np.nan , A_tf)
A_tf_filtered = A_tf_filtered.astype(float)
A_tf_filtered[np.all(np.isnan(A_tf_filtered), axis=1)] = 0

# Compute Mapper
#print(n_tf, n_targets_tf[:,None].shape, np.arange(n_tf)[:,None].shape)
mapper = MapperComplex(input_type='point cloud', min_points_per_node=0, filter_bnds=None, resolutions=[5,5], gains=[0.3,0.3], clustering=grnCC_tf)
# mapper.fit(X=np.arange(n_tf)[:,None], filters=eigenvectors[:,0:2], colors=np.hstack([n_targets_tf[:,None], A_tf[:,1:], F1_tf[:,1:], F2_tf[:,1:]])) #n_targets_tf[:,None])

mapper.fit(X=np.arange(n_tf)[:,None], filters=eigenvectors[:,0:2], colors=np.hstack([n_targets_tf[:,None], A_tf_filtered[:,1:]]))

# Interactive visu
# mapper.save_to_html(file_name="TFs", data_name="GRN", cover_name="uniform", color_name="F1:1")

# 2D visu
mapper_graph = mapper.get_networkx()

print([(v, le.inverse_transform(  [int(i) for i in le_tf.inverse_transform(mapper.node_info_[v]["indices"])]  )) for v in mapper_graph.nodes()])
with open("mapper_nodes.txt", "w") as f:
    result = [(v, le.inverse_transform([int(i) for i in le_tf.inverse_transform(mapper.node_info_[v]["indices"])])) for v in mapper_graph.nodes()]
    f.write(str(result))

mapper_graph.remove_nodes_from(list(nx.isolates(mapper_graph)))

## Custom cmap for tf activity
cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)
# cmap = sns.color_palette('coolwarm', as_cmap=True)
for cond in range(1, A_tf_filtered.shape[1]):
    plt.figure()
    node_colors = [mapper.node_info_[v]["colors"][cond] for v in mapper_graph.nodes()]
    # Specify node sizes and edge weights
    # One possible layout, use the one you prefer
    pos = nx.kamada_kawai_layout(mapper_graph)
    # pos = graphviz_layout(mapper_graph, prog="dot")
    # Specify node colors
    pathcollection = nx.draw_networkx_nodes(mapper_graph, pos, node_color=node_colors, cmap=cmap, vmin=-2.5, vmax=2.5)
    # Specify node sizes and edge weights
    nx.draw(mapper_graph, pos=pos, with_labels=True,
                          node_color=node_colors,
                          cmap=cmap,
                          vmin=-2.5,
                          vmax=2.5,
                          node_size=[10*mapper.node_info_[v]["colors"][0] for v in mapper_graph.nodes()],
                          width=[  common_targets_tf[np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]])),:]
                                                    [:,np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]]))].sum()/100 for e in mapper_graph.edges()])
    # Plot colorscale
    plt.colorbar(pathcollection)
    plt.savefig('./Plot/mapper_filtered_' + str(cond) + '.svg', format='svg', dpi=300, bbox_inches='tight')   

### Plot Graph without coloring
plt.figure()
pos = nx.kamada_kawai_layout(mapper_graph)
nx.draw(mapper_graph, pos=pos, with_labels=False,
                        node_size=[10*mapper.node_info_[v]["colors"][0] for v in mapper_graph.nodes()],
                        width=[  common_targets_tf[np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]])),:]
                                                    [:,np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]]))].sum()/100 for e in mapper_graph.edges()])
### Plot Graph colored by number of TF in Cluster

# signif_tf_cluster = pd.read_csv('./Significant_TFs_per_cluster.tsv', sep='\t', header=0, index_col=0)
# signif_tf_cluster.drop(index=[10,12],inplace=True)
# signif_tf_cluster['Significant'].to_list()


# plt.figure()
# pos = nx.kamada_kawai_layout(mapper_graph)
# number_tf = [len(mapper.node_info_[v]["indices"]) for v in mapper_graph.nodes()]
# pathcollection = nx.draw_networkx_nodes(mapper_graph, pos, node_color=number_tf, cmap='Reds')
# nx.draw(mapper_graph, pos=pos, with_labels=True,
#                         cmap='Reds',
#                         node_color = number_tf,
#                         node_size=[10*mapper.node_info_[v]["colors"][0] for v in mapper_graph.nodes()],
#                         width=[  common_targets_tf[np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]])),:]
#                                                     [:,np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]]))].sum()/100 for e in mapper_graph.edges()])
# plt.colorbar(pathcollection)

# edge_sums = {
#     (e[0], e[1]): common_targets_tf[np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]])), :]
#     [:, np.unique(np.concat([mapper.node_info_[e[0]]["indices"], mapper.node_info_[e[1]]["indices"]]))].sum()
#     for e in mapper_graph.edges()
# }

# # Assuming the number of nodes is known and is `num_nodes`
# num_nodes = len(mapper.node_info_)
# adj_matrix = np.zeros((num_nodes, num_nodes))

# for (i, j), sum_value in edge_sums.items():
#     adj_matrix[i, j] = sum_value
#     adj_matrix[j, i] = sum_value  # Assuming the graph is undirected

# adj_matrix_df = pd.DataFrame(adj_matrix)
# adj_matrix_df_no_zero = adj_matrix_df.loc[~(adj_matrix_df==0).all(axis=1)]
# adj_matrix_df_filtered = adj_matrix_df_no_zero.loc[:, ~(adj_matrix_df_no_zero == 0).all(axis=0)]

# plt.figure(figsize=(12, 10))  # Adjust the figure size for better visibility
# sns.heatmap(
#     adj_matrix,
#     cmap="magma_r",  # Choose a colormap
#     annot=True,      # Enable annotations
#     cbar=True,        # Show the colorbar
#     cbar_kws={'label': 'Number of targets between clusters'},
#     fmt ='g'
# )
# sns.clustermap(
#     adj_matrix_df_filtered,
#     cmap="magma_r",  # Choose a colormap
#     annot=True,      # Enable annotations
#     cbar=True,        # Show the colorbar
#     cbar_kws={'label': 'Number of targets between clusters'},
#     fmt ='g',
#     method='average',
#     metric='correlation'
# )

# adj_matrix_df.to_csv('./adj_matrix.adj', sep='\t', header=True, index=True)
# le.inverse_transform( [int(i) for i in le_tf.inverse_transform([70])])

# #plt.figure()
# #plt.hist([mapper.node_info_[v]["colors"][0] for v in mapper_graph.nodes()], bins=20)
# #plt.show()

# le.inverse_transform( [int(i) for i in le_tf.inverse_transform([70])])
# le_tf.transform( [int(i) for i in le.transform(['Solyc12g099370'])])
# ##split the common targets into cluster
# [mapper.node_info_[v]["indices"] for v in mapper_graph.nodes()]

# ## plot by Clusters
# for v in mapper_graph.nodes():
#     indices = np.array(mapper.node_info_[v]["indices"]) # indices of TF in the cluster node
#     common_targets_tf_save_nodes = common_targets_tf_save[np.ix_(indices, indices)] ##
#     plt.figure(figsize=(16, 12))  # Increase the figure size for better visibility
#     sns.heatmap(
#         common_targets_tf_save_nodes,
#         cmap="magma_r",  # Choose a colormap
#         annot=False,      # Disable annotations for large arrays
#         cbar=True,        # Show the colorbar
#         cbar_kws={'label': 'Shared Targets'},
#         xticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(indices)]), # Even tho we loose it we can access it with indices and transform back into OLN
#         yticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(indices)])
#     )
#     plt.title(f'Adjacency matrix of targets in cluster {v}')
#     # plt.savefig(f'./Adj_shared_targets/Adjacency_targets_cluster_{v}.png', dpi=300, bbox_inches='tight' )   

# cluster = X[np.isin(X[:,0],le_tf.inverse_transform(mapper.node_info_[7]["indices"]))]
# len(np.unique(np.concatenate((cluster[:,0], cluster[:,1]))))

# ##Plot combining clusters
# custom_cluster = [7, 19]
# # custom_cluster = [18,9,13,15,17]
# saved_indices = list()
# for v in custom_cluster:
#     indices = np.array(mapper.node_info_[v]["indices"])  # get the indices of the TF in the cluster
#     saved_indices.extend(indices)
    
# saved_indices = np.unique(saved_indices)

# common_targets_tf_save_nodes = common_targets_tf_save[np.ix_(saved_indices, saved_indices)] ##get the targets of the TF thanks to the indices
# plt.figure(figsize=(16, 12))  # Increase the figure size for better visibility
# sns.heatmap(
#     common_targets_tf_save_nodes,
#     cmap="magma_r",  # Choose a colormap
#     annot=False,      # Disable annotations for large arrays
#     cbar=True,        # Show the colorbar
#     cbar_kws={'label': 'Shared Targets'},
#     xticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(saved_indices)]), # Even tho we loose it we can access it with indices and transform back into OLN
#     yticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(saved_indices)])
# )
# plt.title(f'Adjacency matrix of targets in cluster {custom_cluster}')

# ### Jaccard index of TF targets of combined clusters
# # Idea is to compute the index based on their targets instead of just numbers.

# common_targets_tf_jaccard = np.zeros([n_tf, n_tf])

# for i in range(n_tf):
#     for j in range(i+1,n_tf):
#         common_targets_tf_jaccard[i,j] = len(np.intersect1d(np.array(targets_tf[i]), np.array(targets_tf[j]))) / len(np.union1d(np.array(targets_tf[i]), np.array(targets_tf[j])))
#         common_targets_tf_jaccard[j,i] = len(np.intersect1d(np.array(targets_tf[i]), np.array(targets_tf[j]))) / len(np.union1d(np.array(targets_tf[i]), np.array(targets_tf[j])))

# ## If we want to divide by the total number of targets
# # mega_cluster = X[np.isin(X[:,0],le_tf.inverse_transform(saved_indices))]
# # len(np.unique(np.concatenate((mega_cluster[:,0], mega_cluster[:,1]))))

# common_targets_tf_jaccard_nodes = common_targets_tf_jaccard[np.ix_(saved_indices, saved_indices)]
# plt.figure(figsize=(16, 12))  # Increase the figure size for better visibility
# sns.heatmap(
#     common_targets_tf_jaccard_nodes,
#     cmap="magma_r",  # Choose a colormap
#     annot=False,      # Disable annotations for large arrays
#     cbar=True,        # Show the colorbar
#     cbar_kws={'label': 'Shared Targets'},
#     vmin=0,
#     vmax=0.35,
#     xticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(saved_indices)]), # Even tho we loose it we can access it with indices and transform back into OLN
#     yticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(saved_indices)])
# )
# plt.title(f'Adjacency matrix of targets in cluster {custom_cluster}')

### Function
def jaccard_index_tf_targets(targets_tf, n_tf):
    """
    Calculate the Jaccard index for pairs of transcription factors (TFs) based on their target genes.

    The Jaccard index is a measure of similarity between two sets, defined as the size of the intersection 
    divided by the size of the union of the sets.

    Parameters:
    targets_tf (list of lists): A list where each element is a list of target genes for a specific TF.
    n_tf (int): The number of transcription factors.

    Returns:
    numpy.ndarray: A 2D array of shape (n_tf, n_tf) containing the Jaccard index for each pair of TFs.
    """
    common_targets_tf_jaccard = np.zeros([n_tf, n_tf])

    for i in range(n_tf):
        for j in range(i+1,n_tf):
            common_targets_tf_jaccard[i,j] = len(np.intersect1d(np.array(targets_tf[i]), np.array(targets_tf[j]))) / len(np.union1d(np.array(targets_tf[i]), np.array(targets_tf[j])))
            common_targets_tf_jaccard[j,i] = len(np.intersect1d(np.array(targets_tf[i]), np.array(targets_tf[j]))) / len(np.union1d(np.array(targets_tf[i]), np.array(targets_tf[j])))
    return common_targets_tf_jaccard

def plot_adjacency_matrix(targets_tf, n_tf, custom_cluster):
    common_targets_tf_jaccard = jaccard_index_tf_targets(targets_tf, n_tf)#Compute the Jaccard index of TF targets
    # Get the indices of the TFs in the cluster(s)
    saved_indices = list()
    for v in custom_cluster:
        indices = np.array(mapper.node_info_[v]["indices"])  # get the indices of the TF in the cluster
        saved_indices.extend(indices)
        
    saved_indices = np.unique(saved_indices)
    
    #Clear non significant TFs
    to_remove = list()
    for TF_indices in saved_indices:
        if np.isnan(A_tf_filtered[TF_indices,1:8]).all():
            print(f"TF {le.inverse_transform( [int(i) for i in le_tf.inverse_transform([TF_indices])])} has no significant activity data")
            to_remove = np.append(to_remove, TF_indices)
    saved_indices = np.setdiff1d(saved_indices, to_remove)
    
    # Compute the adjacency matrix of the shared targets based on the global matrix
    common_targets_tf_jaccard_nodes = common_targets_tf_jaccard[np.ix_(saved_indices, saved_indices)]
    
    plt.figure(figsize=(16, 12))  # Increase the figure size for better visibility
    sns.heatmap(
        common_targets_tf_jaccard_nodes,
        cmap="magma_r",  # Choose a colormap
        annot=False,      # Disable annotations for large arrays
        cbar=True,        # Show the colorbar
        cbar_kws={'label': 'Jaccard Index of Shared Targets'},
        vmin=0,
        vmax=0.35,
        xticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(saved_indices)]), # Even tho we loose it we can access it with indices and transform back into OLN
        yticklabels=le.inverse_transform( [int(i) for i in le_tf.inverse_transform(saved_indices)])
    )
    plt.title(f'Adjacency matrix of targets in cluster {custom_cluster}')
    plt.savefig(f'./Adj_shared_targets/Adjacency_targets_cluster_{custom_cluster}.png', dpi=300, bbox_inches='tight' )

plot_adjacency_matrix(targets_tf, n_tf, [0,2,3,5])

### Manually checking
# clusterA = X[np.isin(X[:,0], le_tf.inverse_transform(mapper.node_info_[7]["indices"]))]
# clusterB = X[np.isin(X[:,0], le_tf.inverse_transform(mapper.node_info_[19]["indices"]))]
# combined_cluster = np.concatenate((clusterA, clusterB), axis=0)
# len(np.unique(np.concatenate((combined_cluster[:,0], combined_cluster[:,1]))))


#### Full GRN
#
#n_all = len(le.classes_)
#f1_all,_,all_indices = np.intersect1d(np.arange(n_all), F1[:,0], return_indices=True)
#F1_all = F1[all_indices,:]
#f2_all,_,all_indices = np.intersect1d(np.arange(n_all), F2[:,0], return_indices=True)
#F2_all = F2[all_indices,:]
#
#graph_all = np.zeros([n_all, n_all])
#for pair in X:
#    graph_all[pair[0], pair[1]] = 1
#    graph_all[pair[1], pair[0]] = 1
#    
#grnCC_all = GraphCC(graph_all)
#
#L = nx.normalized_laplacian_matrix(nx.Graph(graph_all))
#eigenvalues, eigenvectors = np.linalg.eig(L.toarray())
#eigenvectors = np.real(eigenvectors)
#print(eigenvectors.shape)
#
#plt.figure()
#plt.hist(eigenvectors[:,0], bins=10)
#plt.show()
#
#mapper = MapperComplex(input_type='point cloud', min_points_per_node=0, filter_bnds=None, resolutions=[5], gains=[0.3], clustering=grnCC_all)
#mapper.fit(X=np.arange(n_all)[:,None], filters=eigenvectors[:,0:1], colors=np.hstack([eigenvectors[:,0:1], F1_all[:,1:], F2_all[:,1:]])) #n_targets_tf[:,None])
#
#mapper.save_to_html(file_name="All", data_name="GRN", cover_name="uniform", color_name="F1:1")
#
#mapper_graph = mapper.get_networkx()
#for cond in range(1,F1_all.shape[1]):
#    plt.figure()
#    nx.draw(mapper_graph, pos=nx.kamada_kawai_layout(mapper_graph), 
#                          node_color=[mapper.node_info_[v]["colors"][cond] for v in mapper_graph.nodes()],
#                          node_size=[mapper.node_info_[v]["colors"][0] for v in mapper_graph.nodes()])
#    plt.show()

