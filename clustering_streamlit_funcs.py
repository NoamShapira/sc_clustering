import collections
import logging
from pathlib import Path
from typing import Tuple, Dict, Callable

import anndata as an
import pandas as pd
import scanpy as sc
import streamlit as st
from matplotlib import pyplot as plt
from sklearn import metrics

from clustering import scanpy_cluster
from clustering.meta_cell import load_meta_cell_and_merge_to_adata
from data import preprocces
from data.data_loading import SeranoDataLoaderFactory, SeranoDataLoaderDescription


def scatter_n_genes_and_n_mt_genes_per_cell(adata, ax_1, ax_2):
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax_1, show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax_2, show=False)


def plot_raw_data(adata):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    sc.pl.highest_expr_genes(adata, n_top=15, ax=ax1, show=False)
    scatter_n_genes_and_n_mt_genes_per_cell(adata, ax2, ax3)
    st.write(fig)


# data loading
@st.cache(allow_output_mutation=True)
def load_data() -> Tuple[an.AnnData, Path]:
    adata, experiment_results_dir_path = SeranoDataLoaderFactory.create_serano_dataloader(
        SeranoDataLoaderDescription.ARM_1_FROM_WEINER).load_data_to_anndata_and_save_to_dir()
    scanpy_cluster.pp_rename_vars_add_mt_metrics(adata)
    return adata, experiment_results_dir_path


# pre-process
@st.cache(allow_output_mutation=True, suppress_st_warning=True)
def drop_bad_genes(adata, filter_labs_bad_genes, min_num_cells_per_gene, min_num_genes_per_cell, max_num_genes_per_cell,
                   max_num_mt_genes_pre_cell, min_num_counts_per_cell):
    st.write("drop_bad_genes cache missed - calculating ... ")
    adata = adata.copy()
    n_genes_source = len(adata.var_names)
    n_cells_source = len(adata.obs_names)
    adata = scanpy_cluster.pp_drop_genes_and_cells(adata,
                                                   min_num_genes_per_cell=min_num_genes_per_cell,
                                                   min_num_cells_per_gene=min_num_cells_per_gene,
                                                   max_num_genes_per_cell=max_num_genes_per_cell,
                                                   max_num_mt_genes_pre_cell=max_num_mt_genes_pre_cell,
                                                   min_num_counts_per_cell=min_num_counts_per_cell)
    n_genes_after_scanpy = len(adata.var_names)
    n_cells_after_scanpy = len(adata.obs_names)
    logging.info(f"scanpy dropped : {n_genes_source - n_genes_after_scanpy} genes, from {n_genes_source}")
    logging.info(f"scanpy dropped : {n_cells_source - n_cells_after_scanpy} cells, from {n_cells_source}")
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
def normalize_and_choose_genes(adata, drop_unvariable_genes, regress_out_total_cont_and_mt):
    st.write("normalize_and_choose_genes cache missed - calculating ... ")
    new_adata = adata.copy()
    scanpy_cluster.normalize_and_choose_genes(new_adata, drop_unvariable_genes, regress_out_total_cont_and_mt)
    return new_adata


# data modeling
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
                                               reference_path=None):
    st.write("computer_clusters_cache cache missed - calculating")
    new_adata = adata.copy()
    scanpy_cluster.cluster_adata(new_adata, method=clustering_method, resolution=resolution)
    if reference_path is not None:
        new_adata = load_meta_cell_and_merge_to_adata(new_adata, reference_path)
    return new_adata


# evaluation
def compute_metrics(adata: an.AnnData, metric_funcs_dict: Dict[str, Callable] = None,
                    clustring_method_name: str = "leiden") -> pd.DataFrame:
    metric_funcs = metric_funcs_dict if metric_funcs_dict is not None else {
        "AMI": metrics.adjusted_mutual_info_score,
        "completness": metrics.completeness_score,
        "ARI": metrics.adjusted_rand_score,
        "homogeneity": metrics.homogeneity_score
    }
    mc_obs = adata.obs["mc.mc"]
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
