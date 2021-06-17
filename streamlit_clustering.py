import collections
import logging
from pathlib import Path
from typing import Tuple, Callable, Dict

import anndata as an
import numpy as np
import pandas as pd
import scanpy as sc
import streamlit as st
from matplotlib import pyplot as plt
from sklearn import metrics

import config
from clustering import scanpy_cluster
from clustering.meta_cell import load_meta_cell_and_merge_to_adata
from data import preprocces
from streamlit_funcs import load_data_and_save_to_results_dir


def compute_metrics(adata: an.AnnData, metric_funcs_dict: Dict[str, Callable] = None,
                    clustring_method_name: str = "leiden") -> pd.DataFrame:
    metric_funcs = metric_funcs_dict if metric_funcs_dict is not None else {
        "AMI": metrics.adjusted_mutual_info_score,
        "completness": metrics.completeness_score,
        "ARI": metrics.adjusted_rand_score,
        "homogeneity": metrics.homogeneity_score
    }
    mc_obs = adata.obs["mc"]
    group_labels = adata.obs["group"]
    sub_group_labels = adata.obs["sub_group"]
    clustering_results = adata.obs[clustring_method_name]
    comperissons = {
        "y_true-mc__y_pred-clustering": (mc_obs, clustering_results),
        "y_true-mc__y_pred-group": (mc_obs, group_labels),
        "y_true-clustering__y_pred-group": (clustering_results, group_labels)
    }
    data = collections.defaultdict(list)
    ind = []
    for comperissons_name, (y_true, y_pred) in comperissons.items():
        ind.append(comperissons_name)
        for metric_name, metric_func in metric_funcs.items():
            data[metric_name].append(metric_func(y_true, y_pred))
    return pd.DataFrame(data=dict(data), index=ind)


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


@st.cache(allow_output_mutation=True, suppress_st_warning=True)
def computer_clusters_cache_and_load_reference(adata, clustering_method, resolution,
                                               reference_path=None, reference_col_name="mc"):
    st.write("computer_clusters_cache cache missed - calculating")
    new_adata = scanpy_cluster.cluster_adata(adata, method=clustering_method, resolution=resolution)
    if reference_path is not None:
        new_adata = load_meta_cell_and_merge_to_adata(new_adata, reference_path, reference_col_name=reference_col_name)
    return new_adata


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

# st.write(sc.pl.pca(adata_4, show=False, return_fig=True))
# sc.pl.pca_variance_ratio(adata, log=True)

st.subheader("build a graph")
n_neighbors = st.number_input("neighborhood graph n neighbors", value=config.NEIGHBORHOOD_GRAPH_N_NEIGHBORS)
n_pcs = st.number_input("neighborhood graph n pcs", value=config.NEIGHBORHOOD_GRAPH_N_PCS)

adata_graph = compute_neighborhood_graph_cache(adata_pca, n_neighbors, n_pcs)
st.write(sc.pl.umap(adata_graph, show=False, return_fig=True))

st.subheader("cluster")
res = st.number_input("clustering resolution", value=config.TL_LEIDEN_RESOLUTION)
clustering_method_name = st.selectbox("Select clustering method", ["leiden"])

reference_col_name = "mc"
final_adata = computer_clusters_cache_and_load_reference(adata_graph, clustering_method=clustering_method_name,
                                                         resolution=res, reference_path=config.META_CELL_PATH,
                                                         reference_col_name=reference_col_name)
st.write(sc.pl.umap(final_adata, ncols=2, show=False, return_fig=True,
                    color=[reference_col_name, clustering_method_name, "group", "sub_group"]))


partition_to_visualize = st.selectbox("choose a partirion of the data to visialize by", [clustering_method_name, reference_col_name, "sub_group"])
chosen_name = st.selectbox(f"select {partition_to_visualize} to show in graph",
                            sorted(final_adata.obs[partition_to_visualize].unique()))
final_adata.obs[f"is_not_chosen_{partition_to_visualize}"] = \
    final_adata.obs[partition_to_visualize].apply(lambda x: x != chosen_name).astype('category')
st.write(sc.pl.umap(final_adata, ncols=1, show=False, return_fig=True,
                    color=f"is_not_chosen_{partition_to_visualize}",
                    palette="Set1"))


st.subheader("Metrict of similarities between partitions")
st.write(compute_metrics(final_adata))

if st.button("save final result to file"):
    final_adata.write(Path(experiment_results_dir_path, "final_adata.h5ad"))
