from pathlib import Path
from typing import Tuple

import anndata as an
import streamlit as st
import scanpy as sc
from matplotlib import pyplot as plt

import config
from clustering import scanpy_cluster
from streamlit_funcs import load_data_and_save_to_results_dir


def scatter_n_genes_and_n_mt_genes_per_cell(ax_1, ax_2):
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax_1, show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax_2, show=False)


@st.cache(allow_output_mutation=True)
def load_data() -> Tuple[an.AnnData, Path]:
    return load_data_and_save_to_results_dir()


adata, experiment_results_dir_path = load_data()

st.write(f"Load data result dir is {experiment_results_dir_path}")
st.title(f"Raw Data")
scanpy_cluster.pp_rename_vars_add_mt_metrics(adata)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
sc.pl.highest_expr_genes(adata, n_top=15, ax=ax1, show=False)
scatter_n_genes_and_n_mt_genes_per_cell(ax2, ax3)
st.write(fig)

st.title("Pre Processing")
st.subheader("drop bad genes and cells")
min_num_cells_per_gene = st.number_input("min number of cells per gene", value=config.PP_MIN_NUMBER_OF_CELLS_PER_GENE)
min_num_genes_per_cell = st.number_input("min number of genes per cell", value=config.PP_MIN_NUMBER_OF_GENES_PER_CELL)
max_num_genes_per_cell = st.number_input("max number of genes per cell", value=config.PP_MAX_NUMBER_OF_GENES_PER_CELL)
max_num_mt_genes_pre_cell = st.number_input("max number of mt genes per cell",
                                            value=config.PP_MAX_NUMBER_OF_MT_GENES_PER_CELL)

adata = scanpy_cluster.pp_drop_genes_and_cells(adata,
                                               min_num_genes_per_cell=min_num_genes_per_cell,
                                               min_num_cells_per_gene=min_num_cells_per_gene,
                                               max_num_genes_per_cell=max_num_genes_per_cell,
                                               max_num_mt_genes_pre_cell=max_num_mt_genes_pre_cell)
fig, (ax1, ax2) = plt.subplots(1, 2)
scatter_n_genes_and_n_mt_genes_per_cell(ax1, ax2)
st.write(fig)

# st.subheader("drop mt or ribo genes")
# drop mt genes or

st.subheader("normalize genes")
scanpy_cluster.normelize_data(adata)  # can add arg - normelized_reads_per_cell

st.subheader("compute variable genes")
scanpy_cluster.compute_variable_genes(adata)  # can add args - "min_mean": 0.0125,"max_mean": 3,"min_disp": 0.5
# fig = plt.Figure()
# fig.add_axes(sc.pl.highly_variable_genes(adata, show=False))
# st.write(fig)

adata = scanpy_cluster.choose_variable_genes(adata)
scanpy_cluster.regress_out_and_scale(adata)

st.subheader("pca")
scanpy_cluster.transform_pca_adata(adata)  # can add - pca_svd_solver, pca_n_comps
st.write(sc.pl.pca(adata, show=False, return_fig=True))
# sc.pl.pca_variance_ratio(adata, log=True)

st.subheader("build a graph")
n_neighbors = st.number_input("neighborhood graph n neighbors", value=config.NEIGHBORHOOD_GRAPH_N_NEIGHBORS)
n_pcs = st.number_input("neighborhood graph n neighbors", value=config.NEIGHBORHOOD_GRAPH_N_PCS)
scanpy_cluster.compute_neighborhood_graph(adata, neighborhood_graph_n_neighbors=n_neighbors,
                                          neighborhood_graph_n_pcs=n_pcs)
st.write(sc.pl.umap(adata, show=False, return_fig=True))

st.subheader("cluster")
res = st.number_input("clustering resolution", value=config.TL_LEIDEN_RESOLUTION)
scanpy_cluster.cluster_adata(adata, method="leiden", resolution=res)
st.write(sc.pl.umap(adata, show=False, return_fig=True, color=['leiden']))
