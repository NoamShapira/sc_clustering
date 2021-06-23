from datetime import datetime
from glob import glob
from pathlib import Path
from typing import Tuple

import numpy as np
from os import listdir
import scanpy as sc
import anndata as ad
from tqdm import tqdm

import config
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


def get_mouse_ingest_evaluation(experiment_name: str = "last", clustering_method: str = "leiden") \
        -> Tuple[dict, dict, dict]:
    clustered_adata = get_last_result_from_experiment(experiment_name=experiment_name)
    clustered_adata.obsp["distances"] = clustered_adata.obsp["distances"].todense()

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


if __name__ == '__main__':
    _, _, _ = get_mouse_ingest_evaluation()
