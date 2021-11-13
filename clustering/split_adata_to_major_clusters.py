import logging
from typing import List

import anndata as ad
import pandas as pd
import scanpy as sc
from sklearn.cluster import MiniBatchKMeans, SpectralClustering, AgglomerativeClustering

import config


def split_to_major_cluster(adata_with_scvi_emb: ad.AnnData, emb_col_name: str,
                           clustering_model_name: str, num_clusters: int) -> List[ad.AnnData]:
    adata_with_scvi_emb_copy = adata_with_scvi_emb.copy()
    col_name_to_split_by = "major_clusters"
    logging.info(f"spliting adata to {num_clusters} clusters, using {clustering_model_name}")
    if clustering_model_name == "kmeans":
        clustering_model = MiniBatchKMeans(n_clusters=num_clusters, verbose=1)
        adata_with_scvi_emb_copy.obs[col_name_to_split_by] = clustering_model.fit(
            adata_with_scvi_emb_copy.obsm[emb_col_name]).labels_
    elif clustering_model_name == "spectral":
        clustering_model = SpectralClustering(n_clusters=num_clusters, affinity='nearest_neighbors', n_jobs=10)
        adata_with_scvi_emb_copy.obs[col_name_to_split_by] = clustering_model.fit(
            adata_with_scvi_emb_copy.obsm[emb_col_name]).labels_
    elif clustering_model_name == "agglomerative":
        clustering_model = AgglomerativeClustering(n_clusters=num_clusters, linkage='ward')
        adata_with_scvi_emb_copy.obs[col_name_to_split_by] = clustering_model.fit(
            adata_with_scvi_emb_copy.obsm[emb_col_name]).labels_
    elif clustering_model_name == "leiden_and_agglomerative":
        clustering_model = AgglomerativeClustering(n_clusters=num_clusters, linkage='ward')

        sc.pp.neighbors(adata_with_scvi_emb_copy, use_rep=emb_col_name)
        leiden_col_name = "leiden"
        sc.tl.leiden(adata_with_scvi_emb_copy, resolution=config.TL_LEIDEN_RESOLUTION * 1.5, key_added=leiden_col_name)

        all_cells_df = pd.DataFrame(index=adata_with_scvi_emb_copy.obs_names,
                                    data=adata_with_scvi_emb_copy.obsm[emb_col_name])
        all_cells_df[leiden_col_name] = adata_with_scvi_emb_copy.obs[leiden_col_name]
        df_of_leiden_means = all_cells_df.groupby(by=leiden_col_name).mean()
        leiden_labels_col_name = "leiden_labels"
        df_of_leiden_means[leiden_labels_col_name] = clustering_model.fit(df_of_leiden_means).labels_
        cells_cluster = pd.merge(left=all_cells_df, right=df_of_leiden_means, how="left",
                                 left_on=leiden_col_name, right_index=True)[leiden_labels_col_name]
        adata_with_scvi_emb_copy.obs[col_name_to_split_by] = cells_cluster
    else:
        raise NotImplementedError

    adatas = []
    for group, idx in adata_with_scvi_emb_copy.obs.groupby(col_name_to_split_by).indices.items():
        sub_adata = adata_with_scvi_emb_copy[idx].copy()
        adatas.append(sub_adata)

    return adatas