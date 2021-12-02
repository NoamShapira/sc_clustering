import anndata as ad
import numpy as np
import scanpy as sc
import streamlit as st
import imageio

import config
from clustering.meta_cell_csv import run_full_pipeline_and_load_meta_cell, MetaCellResultsColumnsNames
from clustering.post_procces_clusters_by_mice import compute_confusion_matrix
from config import CLUSTERING_METHOD
from data.amp_batch_data_loader_factory import PlatesLoaderDescription


@st.cache(allow_output_mutation=True, suppress_st_warning=True)
def full_pipeline_and_load_meta_cell_cached(
        data_loader_desc: PlatesLoaderDescription = PlatesLoaderDescription.ARMS_1_2_3_FROM_NOAMSH,
        path_to_meta_cell_results=config.META_CELL_PATH):
    return run_full_pipeline_and_load_meta_cell(data_loader_desc, path_to_meta_cell_results)

def compute_and_plot_rank_genes_groups_to_local_path(adata, groupby, groups, method, n_genes, prefix,
                                                     *kwargs) -> str:
    sc.tl.rank_genes_groups(adata, groupby=groupby, groups=groups, method=method)
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, save=prefix)
    return f"figures/rank_genes_groups_{groupby}{prefix}"

def visualize_meta_cell(adata, meta_cell):
    tmp_mc_col_name = "meta_cell_to_visualize"
    adata.obs[tmp_mc_col_name] = adata.obs[MetaCellResultsColumnsNames().meta_cell].apply(
        lambda x: 1 if x == meta_cell else 0).astype('category')
    st.write(sc.pl.umap(adata, color=[tmp_mc_col_name], show=False, return_fig=True))
    adata.obs[MetaCellResultsColumnsNames().meta_cell] = adata.obs[MetaCellResultsColumnsNames().meta_cell].astype(
        'str')
    image_path = compute_and_plot_rank_genes_groups_to_local_path(adata,
                                                                  groupby=MetaCellResultsColumnsNames().meta_cell
                                                                  , groups=[meta_cell], method='wilcoxon',
                                                                  n_genes=20, prefix="_fig_1.png")
    st.image(np.asarray(imageio.imread(image_path)))

def get_most_likely_cluster(adata: ad.AnnData, meta_cell_name: str) -> str:
    meta_cells = adata.obs[MetaCellResultsColumnsNames().meta_cell]
    clustering = adata.obs[CLUSTERING_METHOD]
    confusion_matrix = compute_confusion_matrix(meta_cells, clustering)
    return confusion_matrix[meta_cell_name].idxmax()

def visualize_clusters_assosiated_with_meta_cell(adata, meta_cell):
    tmp_cluster_col_name = "cluster_col_name_to_visualize"
    matching_cluster = get_most_likely_cluster(adata, meta_cell)
    adata.obs[tmp_cluster_col_name] = adata.obs[CLUSTERING_METHOD].apply(
        lambda x: 1 if x == matching_cluster else 0).astype('category')
    st.write(sc.pl.umap(adata, color=[tmp_cluster_col_name], show=False, return_fig=True))
    col1, col2 = st.beta_columns(2)
    image_path_1 = compute_and_plot_rank_genes_groups_to_local_path(adata, groupby=CLUSTERING_METHOD,
                                                                    groups=[matching_cluster], method='wilcoxon',
                                                                    n_genes=20, prefix="_fig_2.png")
    col1.header("scanpy cluster vs all")
    col1.image(np.asarray(imageio.imread(image_path_1)))

    adata_only_cluster = adata[adata.obs[CLUSTERING_METHOD] == matching_cluster, :]
    image_path_2 = compute_and_plot_rank_genes_groups_to_local_path(adata_only_cluster,
                                                                    groupby=MetaCellResultsColumnsNames().meta_cell,
                                                                    groups=[meta_cell], method='wilcoxon',
                                                                    n_genes=20, prefix="_fig_3.png")
    col2.header("meta cell vs rest of the scanpy cluster")
    col2.image(np.asarray(imageio.imread(image_path_2)))
