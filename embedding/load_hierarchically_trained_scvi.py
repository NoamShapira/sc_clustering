import os
from pathlib import Path
from typing import Tuple, Dict

import anndata as ad

import config
from scripts.train_scvi_hierarchically import emb_col_name


def load_adatas_from_experiment(experiment_name: str, result_dir: Path = config.RESULTS_DIR) -> \
        Tuple[ad.AnnData, ad.AnnData, Dict[str, ad.AnnData], Dict[str, str]]:
    experiment_dir = Path(result_dir, experiment_name)
    full_adata_file_name = f"adata_with_embedding_in_{emb_col_name}.h5ad"
    adata_with_scanpy = ad.read_h5ad(Path(experiment_dir, "scanpy_metacells.h5ad"))
    full_adata = ad.read_h5ad(Path(experiment_dir, full_adata_file_name))
    clusters_adatas = {}
    emb_coll_names = {}
    for root, subdirectories, _ in os.walk(experiment_dir):
        for sub_dir in [s_dir for s_dir in subdirectories if "cluster" in s_dir]:
            for file in [f for f in os.listdir(Path(root, sub_dir)) if ".h5ad" in f]:
                clusters_adatas[sub_dir] = ad.read_h5ad(Path(root, sub_dir, file))
                emb_coll_names[sub_dir] = file.split(".")[0].split("embedding_in_")[1]
    return full_adata, adata_with_scanpy, clusters_adatas, emb_coll_names

# full_adata, adata_with_scanpy, clusters_adatas, emb_coll_names = load_adatas_from_experiment('hierarchical_scvi_2021_11_18__10_10_08')
