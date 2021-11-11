import argparse
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering, SpectralClustering, MiniBatchKMeans

sys.path.append('/home/labs/amit/noamsh/repos/sc_clustering')
import config
from data.meta_data_columns_names import TREATMENT_ARM
from embedding.scvi_pipe import train_scvi_on_adata
from utils import get_now_timestemp_as_string


def split_to_major_cluster(adata_with_scvi_emb: ad.AnnData, emb_col_name: str,
                           clustering_model_name: str = "kmeans",
                           num_clusters: int = 2) -> List[ad.AnnData]:
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
        sc.tl.leiden(adata_with_scvi_emb_copy, resolution=config.TL_LEIDEN_RESOLUTION)

        all_cells_df = pd.DataFrame(index=adata_with_scvi_emb_copy.obs_names,
                                    data=adata_with_scvi_emb_copy.obsm[emb_col_name])
        df_of_leiden_means =
        leiden_labels = clustering_model.fit(
            df_of_leiden_means).labels_

        all_cells_df["leiden"] = adata_with_scvi_emb_copy.obs["leiden"]

        adata_with_scvi_emb_copy[col_name_to_split_by] =

    else:
        raise NotImplementedError

    adatas = []
    for group, idx in adata_with_scvi_emb_copy.obs.groupby(col_name_to_split_by).indices.items():
        sub_adata = adata_with_scvi_emb_copy[idx].copy()
        adatas.append(sub_adata)

    return adatas


experiment_name = f"hierarchical_scvi_{get_now_timestemp_as_string()}"
experiment_results_dir_path = Path(config.RESULTS_DIR, experiment_name)
asaf_raw_adata_path = Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "cells.h5ad")
new_adata_path = Path(experiment_results_dir_path, "cells.h5ad")
asaf_labeled_adata_path = Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "scanpy_metacells.h5ad")
new_labeled_adata_path = Path(experiment_results_dir_path, "scanpy_metacells.h5ad")

parser = argparse.ArgumentParser()
parser.add_argument('--test_mod', type=bool, default=False)
parser.add_argument("--load_raw_adata_from", type=str, default=str(asaf_raw_adata_path))
parser.add_argument('--save_raw_adata_to', type=str, default=str(new_adata_path))
parser.add_argument("--load_labeled_adata_from", type=str, default=str(asaf_labeled_adata_path))
parser.add_argument("--save_labeled_adata_to", type=str, default=str(new_labeled_adata_path))

parser.add_argument("--split_and_continue_from_dir", type=str, default="")

parser.add_argument("--pp_n_top_variable_genes", type=int, default=3000)
parser.add_argument("--clustering_method", type=str, default="kmeans")

parser.add_argument("--n_layers", type=int, default=1)
parser.add_argument("--n_latent", type=int, default=30)
parser.add_argument("--dropout_rate", type=float, default=0.5)
parser.add_argument("--n_hidden", type=int, default=64)

args = parser.parse_args()

os.mkdir(experiment_results_dir_path)
shutil.copyfile(args.load_raw_adata_from, args.save_raw_adata_to)
shutil.copyfile(args.load_labeled_adata_from, args.save_labeled_adata_to)

adata = ad.read_h5ad(args.save_raw_adata_to)
# adata_with_labels = ad.read_h5ad(args.save_labeled_adata_to)
if args.test_mod:
    adata = adata[0:300, :]
    # adata_with_labels = adata_with_labels[adata.obs, :]

emb_col_name = "X_scVI"

if args.split_and_continue_from_dir == "":
    adata_with_scvi_emb = train_scvi_on_adata(adata, args_for_scvi_model=args, results_dir=experiment_results_dir_path,
                                              return_adata_with_embedding=True, emb_col_name=emb_col_name,
                                              batch_col_name=TREATMENT_ARM,
                                              n_variable_genes=args.pp_n_top_variable_genes)
else:
    adata_with_scvi_emb = ad.read_h5ad(
        Path(args.split_and_continue_from_dir, f"adata_with_embedding_in_{emb_col_name}.h5ad"))
    adata_with_scvi_emb.write(Path(experiment_results_dir_path, f"loaded_adata_with_embedding_in_{emb_col_name}.h5ad"))

adatas_of_major_cluster = split_to_major_cluster(adata_with_scvi_emb, emb_col_name=emb_col_name,
                                                 clustering_model_name=args.clustering_method)

for i, adata_of_cluster in enumerate(adatas_of_major_cluster):
    split_results_dir_path = Path(experiment_results_dir_path, f"cluster_{i}")
    os.mkdir(split_results_dir_path)
    train_scvi_on_adata(adata_of_cluster, args_for_scvi_model=args,
                        results_dir=split_results_dir_path, return_adata_with_embedding=True,
                        emb_col_name=f"{emb_col_name}_{i}", batch_col_name=TREATMENT_ARM,
                        n_variable_genes=args.pp_n_top_variable_genes)

print("finished script !!")
