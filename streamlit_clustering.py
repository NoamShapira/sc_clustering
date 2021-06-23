from pathlib import Path

import scanpy as sc
import streamlit as st

import config
from streamlit_funcs import compute_metrics, load_data, plot_raw_data, \
    drop_bad_genes, normalize_and_choose_genes, compute_pca, compute_neighborhood_graph_cache, \
    computer_clusters_cache_and_load_reference
from utils import get_now_timestemp_as_string

raw_adata, experiment_results_dir_path = load_data()

st.write(f"Load data result dir is {experiment_results_dir_path}")
st.title(f"Raw Data")
plot_raw_data(raw_adata)

st.title("Pre Processing")
st.subheader("drop bad genes and cells")
min_num_cells_per_gene = st.number_input("min number of cells per gene", value=config.PP_MIN_NUMBER_OF_CELLS_PER_GENE)
min_num_genes_per_cell = st.number_input("min number of genes per cell", value=config.PP_MIN_NUMBER_OF_GENES_PER_CELL)
max_num_genes_per_cell = st.number_input("max number of genes per cell", value=config.PP_MAX_NUMBER_OF_GENES_PER_CELL)
min_num_counts_per_cell = st.number_input("min number of counts per cell", value=config.PP_MIN_NUMBER_OF_GENES_PER_CELL)

max_num_mt_genes_pre_cell = st.number_input("max percentage of mt genes per cell",
                                            value=config.PP_MAX_PCT_OF_MT_GENES_PER_CELL)

filter_lab_genes = st.checkbox("Filter bad genes from labs black list")
adata_dropped_genes = drop_bad_genes(raw_adata, filter_lab_genes, min_num_cells_per_gene, min_num_genes_per_cell,
                                     max_num_genes_per_cell, max_num_mt_genes_pre_cell, min_num_counts_per_cell)

st.subheader("normalize and choose genes")
drop_unvariable_genes = st.checkbox("choose variable genes", value=True)
regress_out_total_cont_and_mt = st.checkbox("regress out total gene count and mitochondrial genes", value=True)
adata_norm = normalize_and_choose_genes(adata_dropped_genes, drop_unvariable_genes, regress_out_total_cont_and_mt)

st.subheader("pca")
adata_pca = compute_pca(adata_norm)

st.subheader("build a graph")
n_neighbors = st.number_input("neighborhood graph n neighbors", value=config.NEIGHBORHOOD_GRAPH_N_NEIGHBORS)
n_pcs = st.number_input("neighborhood graph n pcs", value=config.NEIGHBORHOOD_GRAPH_N_PCS)

adata_graph = compute_neighborhood_graph_cache(adata_pca, n_neighbors, n_pcs)
st.write(sc.pl.umap(adata_graph, show=False, return_fig=True))

st.subheader("cluster")
res = st.number_input("clustering resolution", value=config.TL_LEIDEN_RESOLUTION)
clustering_method_name = st.selectbox("Select clustering method", ["leiden"])

reference_col_names = ["mc.mc", "group", "sub_group"]
final_adata = computer_clusters_cache_and_load_reference(adata_graph, clustering_method=clustering_method_name,
                                                         resolution=res, reference_path=config.META_CELL_PATH,
                                                         reference_col_names=reference_col_names)
st.write(sc.pl.umap(final_adata, ncols=2, show=False, return_fig=True,
                    color=[clustering_method_name] + reference_col_names))

partition_to_visualize = st.selectbox("choose a partirion of the data to visialize by",
                                      [clustering_method_name] + reference_col_names)
final_adata_for_visualization = final_adata.copy()
chosen_name = st.selectbox(f"select {partition_to_visualize} to show in graph",
                           sorted(final_adata_for_visualization.obs[partition_to_visualize].unique()))
final_adata_for_visualization.obs[f"is_not_chosen_{partition_to_visualize}"] = \
    final_adata_for_visualization.obs[partition_to_visualize].apply(lambda x: x != chosen_name).astype('category')
st.write(sc.pl.umap(final_adata_for_visualization, ncols=1, show=False, return_fig=True,
                    color=f"is_not_chosen_{partition_to_visualize}",
                    palette="Set1"))

st.subheader("Metrict of similarities between partitions")
st.write(compute_metrics(final_adata))

if st.button("save final result to file"):
    final_adata.write(Path(experiment_results_dir_path, f"final_adata_{get_now_timestemp_as_string()}.h5ad"))
