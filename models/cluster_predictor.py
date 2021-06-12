from collections import Counter
from typing import Tuple, List, Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import spatial

from clustering.scanpy_cluster import normelize_data


class OnlineClustering:
    def __init__(self, adata: ad.AnnData, clustering_obs_col_name: str):
        self.adata = adata
        self.source_adata = adata.copy()
        self.predicted_obs = []
        self.clustering_obs_col_name = clustering_obs_col_name
        self.distance_tree = spatial.KDTree(self.adata.X)

    def drop_predicted_cells(self):
        self.adata = self.source_adata.copy()
        self.predicted_obs = []

    def predict(self, adata_to_pred: ad.AnnData, clusters_to_keep_together_col_in_adata: str = None,
                predictions_col: str = "cluster_predictions"):
        NO_CLUSTER = "NO_CLUSTER"
        if clusters_to_keep_together_col_in_adata is None:
            clusters_to_keep_together_col_in_adata = "clusters_to_assign_to_the_same_group"
            adata_to_pred.obs[clusters_to_keep_together_col_in_adata] = np.arange(len(adata_to_pred.obs))
        else:
            adata_to_pred.obs[clusters_to_keep_together_col_in_adata].fillna(NO_CLUSTER)

        obs_name_to_predictions = []
        for name, obs_group in adata_to_pred.obs.groupby(clusters_to_keep_together_col_in_adata):
            if name == NO_CLUSTER:
                for obs_name in obs_group.index:
                    obs_name_to_predictions.append(self.assign_to_single_cluster(adata_to_pred[obs_name, :]))
            else:
                obs_name_to_predictions.append(self.assign_to_single_cluster(adata_to_pred[obs_group.index, :]))
        obs_name_to_predictions = [item for sublist in obs_name_to_predictions for item in sublist]
        series_of_predictions = pd.Series([x[1] for x in obs_name_to_predictions],
                                          index=[x[0] for x in obs_name_to_predictions])
        adata_to_pred.obs[predictions_col] = series_of_predictions

    def assign_to_single_cluster(self, adata_to_pred: ad.AnnData, method: str = "nearest_neighbor_majority_vote") -> \
            List[Tuple[str, Any]]:
        adata_to_pred = adata_to_pred.copy()
        if method == "nearest_neighbor_majority_vote":
            predicted_clusters = [self.get_nearest_neighbor(adata_to_pred[observation, :])[1] for observation in
                                  adata_to_pred.obs_names]
            predicted_cluster = Counter(predicted_clusters).most_common(1)[0][0]
        elif method == "svm":
            raise NotImplemented
            # predicted_cluster =
        else:
            raise Exception("Ileargal method in - assign_to_single_cluster")

        adata_to_pred.obs["clusters"] = predicted_cluster
        self.predicted_obs.append(adata_to_pred.obs_names)
        return [(obs_name, predicted_cluster) for obs_name in adata_to_pred.obs_names]

    def get_nearest_neighbor(self, adata_to_pred: ad.AnnData, normalize: bool = True) -> Tuple[str, int]:
        query_adata = adata_to_pred.copy()
        if normalize:
            normelize_data(query_adata)
        _, ind_in_adata_X = self.distance_tree.query(query_adata.X.squeeze())
        nearest_neighbor_name = self.adata.obs_names[ind_in_adata_X]
        nearest_neighbor_cluster = self.adata.obs[self.clustering_obs_col_name][nearest_neighbor_name]

        return nearest_neighbor_name, nearest_neighbor_cluster
