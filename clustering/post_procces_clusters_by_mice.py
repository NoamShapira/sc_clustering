import logging
import pickle
from datetime import datetime
from functools import partial
from glob import glob
from os import listdir
from pathlib import Path
from typing import Tuple

import plotly.graph_objects as go
import anndata as ad
import numpy as np
import scanpy as sc
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

import config
from clustering import scanpy_cluster
from data.meta_data_columns_names import ARM_DAY_BATCHMOUSE
from utils import TIME_PATTERN, TIME_PATTERN_LEN


def get_name_of_last_experiment() -> str:
    all_experiment_names = listdir(config.RESULTS_DIR)
    all_experimnt_times = [datetime.strptime(exp_name[-TIME_PATTERN_LEN:], TIME_PATTERN) for exp_name in
                           all_experiment_names]
    return all_experiment_names[np.argmax(all_experimnt_times)]


def get_last_result_from_experiment(experiment_name: str) -> ad.AnnData:
    if experiment_name == "last":
        experiment_name = get_name_of_last_experiment()

    experiment_path = Path(config.RESULTS_DIR, experiment_name)
    all_results_paths = [Path(path) for path in glob(str(experiment_path) + "/final*")]
    all_results_times = [datetime.strptime(result_path.name.split(".")[0][-TIME_PATTERN_LEN:], TIME_PATTERN)
                         for result_path in all_results_paths]
    return ad.read(all_results_paths[np.argmax(all_results_times)])


def get_raw_data_from_experiment(experiment_name: str) -> ad.AnnData:
    if experiment_name == "last":
        experiment_name = get_name_of_last_experiment()

    return ad.read(Path(config.RESULTS_DIR, experiment_name, "loaded_data.h5ad"))


def get_mouse_ingest_evaluation_no_reclustering(experiment_name: str = "last", clustering_method: str = "leiden") \
        -> Tuple[dict, dict, dict]:
    clustered_adata = get_last_result_from_experiment(experiment_name=experiment_name)

    mouse_col_name = "mouse"
    clustered_adata.obs[mouse_col_name] = clustered_adata.obs[ARM_DAY_BATCHMOUSE].apply(lambda x: x.split(".")[-1])
    all_mice = clustered_adata.obs[mouse_col_name].unique()

    evaluation = {}
    clusters_counts_original = {}
    clusters_counts_ingested = {}
    clusters_counts_original["all_mice"] = np.unique(clustered_adata.obs[clustering_method], return_counts=True)

    for mouse in tqdm(all_mice):
        mouse_ind = clustered_adata.obs[mouse_col_name].apply(lambda x: x == mouse)
        ref_adata = clustered_adata[~mouse_ind, :]
        mouse_adata = clustered_adata[mouse_ind, :]

        original_clusters = mouse_adata.obs[clustering_method]
        sc.tl.ingest(adata=mouse_adata, adata_ref=ref_adata, obs=clustering_method)
        ingested_clusters = mouse_adata.obs[clustering_method]
        evaluation[mouse] = (original_clusters, ingested_clusters)
        clusters_counts_original[mouse] = np.unique(original_clusters, return_counts=True)
        clusters_counts_ingested[mouse] = np.unique(ingested_clusters, return_counts=True)

    return evaluation, clusters_counts_original, clusters_counts_ingested


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

    list_of_tuples = process_map(partial(evaluate_ingest_for_mouse, adata=adata, clustered_adata=clustered_adata,
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


def evaluate_ingest_for_mouse(mouse, adata, clustered_adata, clustering_method, mouse_col_name):
    mouse_ind = adata.obs[mouse_col_name].apply(lambda x: x == mouse)
    ref_adata = adata[~mouse_ind, :]
    mouse_adata = adata[mouse_ind, :]
    logging.info(f"clustering ref for mouse {mouse} ...")
    ref_adata = scanpy_cluster.run_full_pipe_from_config(ref_adata, filter_cells_only=False)
    mouse_adata = scanpy_cluster.pp_choose_genes_and_normelize(mouse_adata, filter_cells_only=True)
    mouse_adata = mouse_adata[:, ref_adata.var_names]
    logging.info(f"ingesting for mouse {mouse} ...")
    sc.tl.ingest(adata=mouse_adata, adata_ref=ref_adata, obs=clustering_method)
    ingested_clusters = mouse_adata.obs[clustering_method]
    mouse_cells_in_clustered = [ind for ind in adata[mouse_ind, :].obs_names if ind in clustered_adata.obs_names]
    original_clusters = clustered_adata[mouse_cells_in_clustered, :].obs[clustering_method]
    return mouse, ingested_clusters, original_clusters


def create_sankey_graph_for_clustering(original_clusters, ingested_clusters):
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


if __name__ == '__main__':
    evaluation, clusters_counts_original, clusters_counts_ingested = \
        evaluate_mouse_ingest_error(save_results_as_pickle=True)
