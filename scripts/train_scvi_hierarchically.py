import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import List

import anndata as ad
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering

import config
from data.meta_data_columns_names import TREATMENT_ARM
from embedding.scvi_pipe import train_scvi_on_adata
from utils import get_now_timestemp_as_string

sys.path.append(str(Path(config.AmitLab_Path, 'noamsh')))

def split_to_major_cluster(adata_with_scvi_emb: ad.AnnData, emb_col_name:str,
                           clustering_model_name="agglomerative") -> List[ad.AnnData]:
    adatas = []
    col_name_to_split_by = "major_clusters"
    if clustering_model_name == "agglomerative":
        clustering_model = AgglomerativeClustering(n_clusters=2, linkage='ward')
        adata_with_scvi_emb.obs[col_name_to_split_by] = clustering_model.fit(adata.obsm[emb_col_name]).labels_
    else:
        raise NotImplementedError

    for group, idx in adata_with_scvi_emb.obs.groupby(col_name_to_split_by).indices.items():
        sub_adata = adata_with_scvi_emb[idx].copy()
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
parser.add_argument("--pp_n_top_variable_genes", type=int, default=3000)

args = parser.parse_args()

os.mkdir(experiment_results_dir_path)
shutil.copyfile(args.load_raw_adata_from, args.save_raw_adata_to)
shutil.copyfile(args.load_labeled_adata_from, args.save_labeled_adata_to)

adata = ad.read_h5ad(args.save_raw_adata_to)
adata_with_labels = ad.read_h5ad(args.save_labeled_adata_to)
if args.test_mod:
    adata = adata[0:500, :]
    adata_with_labels = adata_with_labels[adata.obs, :]

adata_with_scvi_emb = train_scvi_on_adata(adata, args_for_scvi_model=args, results_dir=experiment_results_dir_path,
                                          return_adata_with_embedding=True, emb_col_name="X_scVI",
                                          batch_col_name=TREATMENT_ARM, n_variable_genes=args.pp_n_top_variable_genes)
adatas_of_major_cluster = split_to_major_cluster(adata_with_scvi_emb, emb_col_name="X_scVI")

for i, adata_of_cluster in enumerate(adatas_of_major_cluster):
    split_results_dir_path = Path(experiment_results_dir_path, f"cluster_{i}")
    cur_adata_with_scvi_emb = train_scvi_on_adata(adata_of_cluster, args_for_scvi_model=args,
                                                  results_dir=split_results_dir_path, return_adata_with_embedding=True,
                                                  emb_col_name=f"X_scVI_{i}", batch_col_name=TREATMENT_ARM,
                                                  n_variable_genes=args.pp_n_top_variable_genes)
