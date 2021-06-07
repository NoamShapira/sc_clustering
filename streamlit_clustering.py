import logging
from pathlib import Path
from typing import Tuple

import anndata as an
import pandas as pd
import streamlit as st
import scanpy as sc
from sklearn import metrics
from matplotlib import pyplot as plt

import config
from clustering import scanpy_cluster
from clustering.meta_cell import load_meta_cell_and_merge_to_adata
from data import preprocces
from streamlit_funcs import load_data_and_save_to_results_dir


def scatter_n_genes_and_n_mt_genes_per_cell(adata, ax_1, ax_2):
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax_1, show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax_2, show=False)


@st.cache(allow_output_mutation=True)
def load_data() -> Tuple[an.AnnData, Path]:
    adata, experiment_results_dir_path = load_data_and_save_to_results_dir()
    scanpy_cluster.pp_rename_vars_add_mt_metrics(adata)
    return adata, experiment_results_dir_path


def plot_raw_data(adata):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    sc.pl.highest_expr_genes(adata, n_top=15, ax=ax1, show=False)
    scatter_n_genes_and_n_mt_genes_per_cell(adata, ax2, ax3)
    st.write(fig)


@st.cache(allow_output_mutation=True, suppress_st_warning=True)
def drop_bad_genes(adata, filter_labs_bad_genes, min_num_cells_per_gene, min_num_genes_per_cell, max_num_genes_per_cell,
                   max_num_mt_genes_pre_cell):
    st.write("drop_bad_genes cache missed - calculating ... ")
    adata = adata.copy()
    n_genes_source = len(adata.var_names)
    adata = scanpy_cluster.pp_drop_genes_and_cells(adata,
                                                   min_num_genes_per_cell=min_num_genes_per_cell,
                                                   min_num_cells_per_gene=min_num_cells_per_gene,
                                                   max_num_genes_per_cell=max_num_genes_per_cell,
                                                   max_num_mt_genes_pre_cell=max_num_mt_genes_pre_cell)
    n_genes_after_scanpy = len(adata.var_names)
    logging.info(f"scanpy dropped : {n_genes_source - n_genes_after_scanpy} genes, from {n_genes_source}")
    fig, (ax1, ax2) = plt.subplots(1, 2)
    scatter_n_genes_and_n_mt_genes_per_cell(adata, ax1, ax2)
    st.write(fig)
    if filter_labs_bad_genes:
        adata = preprocces.drop_bad_genes(adata)
        n_genes_after_scanpy_and_labs_blacklist = len(adata.var_names)
        logging.info(f"lab blacklist dropped : {n_genes_after_scanpy - n_genes_after_scanpy_and_labs_blacklist}"
                     f" genes, from {n_genes_source}")
    return adata


@st.cache(allow_output_mutation=True, suppress_st_warning=True)
def normalize_and_choose_genes(adata):
    st.write("normalize_and_choose_genes cache missed - calculating ... ")
    return scanpy_cluster.normalize_and_choose_genes(adata)
    # scanpy_cluster.normelize_data(adata)  # can add arg - normelized_reads_per_cell
    # scanpy_cluster.compute_variable_genes(adata)  # can add args - "min_mean": 0.0125,"max_mean": 3,"min_disp": 0.5
    # # fig = plt.Figure()
    # # fig.add_axes(sc.pl.highly_variable_genes(adata, show=False))
    # # st.write(fig)
    #
    # adata = scanpy_cluster.choose_variable_genes(adata)
    # scanpy_cluster.regress_out_and_scale(adata)
    # return adata.copy()


@st.cache(allow_output_mutation=True, suppress_st_warning=True)
def compute_pca(adata):
    st.write("compute_pca cache missed - calculating ... ")
    new_adata = adata.copy()
    scanpy_cluster.transform_pca_adata(new_adata)  # can add - pca_svd_solver, pca_n_comps
    return new_adata


@st.cache(allow_output_mutation=True, suppress_st_warning=True)
def compute_neighborhood_graph_cache(adata, n_neighbors, n_pcs):
    st.write("compute_neighborhood_graph_cache cache missed - calculating")
    copy_adata = adata.copy()
    scanpy_cluster.compute_neighborhood_graph(copy_adata, neighborhood_graph_n_neighbors=n_neighbors,
                                              neighborhood_graph_n_pcs=n_pcs)
    return copy_adata


raw_adata, experiment_results_dir_path = load_data()

st.write(f"Load data result dir is {experiment_results_dir_path}")
st.title(f"Raw Data")
plot_raw_data(raw_adata)

st.title("Pre Processing")
st.subheader("drop bad genes and cells")
min_num_cells_per_gene = st.number_input("min number of cells per gene", value=config.PP_MIN_NUMBER_OF_CELLS_PER_GENE)
min_num_genes_per_cell = st.number_input("min number of genes per cell", value=config.PP_MIN_NUMBER_OF_GENES_PER_CELL)
max_num_genes_per_cell = st.number_input("max number of genes per cell", value=config.PP_MAX_NUMBER_OF_GENES_PER_CELL)

max_num_mt_genes_pre_cell = st.number_input("max number of mt genes per cell",
                                            value=config.PP_MAX_NUMBER_OF_MT_GENES_PER_CELL)

filter_lab_genes = st.checkbox("Filter bad genes from labs black list")
adata_dropped_genes = drop_bad_genes(raw_adata, filter_lab_genes, min_num_cells_per_gene, min_num_genes_per_cell,
                                     max_num_genes_per_cell, max_num_mt_genes_pre_cell)

st.subheader("normalize and choose genes")
adata_norm = normalize_and_choose_genes(adata_dropped_genes)

st.subheader("pca")
adata_pca = compute_pca(adata_norm)

adata_4 = load_meta_cell_and_merge_to_adata(adata_pca, config.META_CELL_PATH)

# st.write(sc.pl.pca(adata_4, show=False, return_fig=True))
# sc.pl.pca_variance_ratio(adata, log=True)

st.subheader("build a graph")
n_neighbors = st.number_input("neighborhood graph n neighbors", value=config.NEIGHBORHOOD_GRAPH_N_NEIGHBORS)
n_pcs = st.number_input("neighborhood graph n pcs", value=config.NEIGHBORHOOD_GRAPH_N_PCS)

adata_graph = compute_neighborhood_graph_cache(adata_4, n_neighbors, n_pcs)
st.write(sc.pl.umap(adata_graph, show=False, return_fig=True))

st.subheader("cluster")
res = st.number_input("clustering resolution", value=config.TL_LEIDEN_RESOLUTION)
final_adata = scanpy_cluster.cluster_adata(adata_graph, method="leiden", resolution=res)
st.write(sc.pl.umap(final_adata, show=False, return_fig=True, color=['leiden', "mc"]))

st.subheader("Metrict of similarities between partitions")
st.write(pd.DataFrame.from_dict({
    "AMI": [metrics.adjusted_mutual_info_score(final_adata.obs["mc"], final_adata.obs["leiden"])],
    "completness": [metrics.completeness_score(final_adata.obs["mc"], final_adata.obs["leiden"])]
}))

if st.button("save final result to file"):
    final_adata.write(Path(experiment_results_dir_path, "final_adata.h5ad"))
