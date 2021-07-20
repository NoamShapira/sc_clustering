import logging
import pickle
from functools import partial
from pathlib import Path
from typing import Tuple

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
import seaborn as sns
from tqdm.contrib.concurrent import process_map

import config
from clustering import scanpy_cluster
from clustering.clustering_experiment_loading import get_name_of_last_experiment, get_last_result_from_experiment, \
    get_raw_data_from_experiment
from data.meta_data_columns_names import ARM_DAY_BATCHMOUSE


def evaluate_mouse_ingest_error(experiment_name: str = "last", clustering_method: str = "leiden",
                                save_results_as_pickle: bool = False):
    adata = get_raw_data_from_experiment(experiment_name)
    clustered_adata = get_last_result_from_experiment(experiment_name)

    mouse_col_name = "mouse"
    adata.obs[mouse_col_name] = adata.obs[ARM_DAY_BATCHMOUSE].apply(lambda x: x.split(".")[-1])
    all_mice = adata.obs[mouse_col_name].unique()

    evaluation = {}
    clusters_counts_original = {}
    clusters_counts_ingested = {}
    clusters_counts_original["all_mice"] = np.unique(clustered_adata.obs[clustering_method], return_counts=True)

    list_of_tuples = process_map(partial(_evaluate_ingest_for_mouse, adata=adata, clustered_adata=clustered_adata,
                                         clustering_method=clustering_method, mouse_col_name=mouse_col_name),
                                 all_mice, unit="mouse", max_workers=10)

    for mouse, ingested_clusters, original_clusters in list_of_tuples:
        evaluation[mouse] = (original_clusters, ingested_clusters)
        clusters_counts_original[mouse] = np.unique(original_clusters, return_counts=True)
        clusters_counts_ingested[mouse] = np.unique(ingested_clusters, return_counts=True)

    if save_results_as_pickle:
        experiment_name = experiment_name if experiment_name != "last" else get_name_of_last_experiment()
        path = Path(config.RESULTS_DIR, experiment_name, "mouse_ingest_results.pkl")
        pickle.dump((evaluation, clusters_counts_original, clusters_counts_ingested), open(path, "wb"))

    return evaluation, clusters_counts_original, clusters_counts_ingested


def _evaluate_ingest_for_mouse(mouse: str, adata: ad.AnnData, clustered_adata: ad.AnnData, clustering_method: str,
                               mouse_col_name: str) -> Tuple[str, pd.Series, pd.Series]:
    mouse_ind = adata.obs[mouse_col_name].apply(lambda x: x == mouse)
    ref_adata = adata[~mouse_ind, :]
    mouse_adata = adata[mouse_ind, :]

    logging.info(f"clustering ref for mouse {mouse} ...")
    ref_adata = scanpy_cluster.run_full_pipe_from_config(ref_adata, filter_cells_only_during_pp=False)
    mouse_adata = scanpy_cluster.pp_choose_genes_and_normelize(mouse_adata, filter_cells_only_during_pp=True)
    mouse_adata = mouse_adata[:, ref_adata.var_names]

    logging.info(f"ingesting for mouse {mouse} ...")
    sc.tl.ingest(adata=mouse_adata, adata_ref=ref_adata, obs=clustering_method)

    ingested_clusters = mouse_adata.obs[clustering_method]
    mouse_cells_in_clustered = [ind for ind in adata[mouse_ind, :].obs_names if ind in clustered_adata.obs_names]
    original_clusters = clustered_adata[mouse_cells_in_clustered, :].obs[clustering_method]
    return mouse, ingested_clusters, original_clusters


def create_sankey_graph_for_clustering(original_clusters: pd.Series, ingested_clusters: pd.Series) -> go.Figure:
    unique_orig_clusters = list(np.unique(original_clusters))
    unique_ingest_clusters = list(np.unique(ingested_clusters))

    links = {"source": [], "target": [], "value": []}
    num_original_clusters = len(unique_orig_clusters)
    for i, orig_cluster in enumerate(unique_orig_clusters):
        for j, ingest_cluster in enumerate(unique_ingest_clusters):
            orig_cells = set(original_clusters[original_clusters == orig_cluster].index)
            ingest_cells = set(ingested_clusters[ingested_clusters == ingest_cluster].index)
            intersecting_cells = orig_cells.intersection(ingest_cells)
            if len(intersecting_cells) > 0:
                links["source"].append(i)
                links["target"].append(j + num_original_clusters)
                links["value"].append(len(intersecting_cells))

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=10,
            line=dict(color="black", width=0.5),
            label=unique_orig_clusters + unique_ingest_clusters,
            color="blue"
        ),
        link=links)])
    return fig


def create_heat_confusion_map(original_clustering: pd.Series, new_clustering: pd.Series):
    confusion_matrix = compute_confusion_matrix(new_clustering, original_clustering)

    return sns.heatmap(confusion_matrix)


def compute_confusion_matrix(new_clustering: pd.Series, original_clustering: pd.Series) -> pd.DataFrame:
    unique_orig_clusters = list(np.unique(original_clustering))
    unique_new_clusters = list(np.unique(new_clustering))
    confusion_matrix = np.zeros((len(unique_orig_clusters), len(unique_new_clusters)))
    for i, orig_cluster in enumerate(unique_orig_clusters):
        for j, new_cluster in enumerate(unique_new_clusters):
            orig_cells = set(original_clustering[original_clustering == orig_cluster].index)
            ingest_cells = set(new_clustering[new_clustering == new_cluster].index)
            confusion_matrix[i, j] = len(orig_cells.intersection(ingest_cells))

    return pd.DataFrame(confusion_matrix, index=unique_orig_clusters, columns=unique_new_clusters)


if __name__ == '__main__':
    _, _, _ = evaluate_mouse_ingest_error(save_results_as_pickle=True)
