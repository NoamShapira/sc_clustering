import os
from pathlib import Path
from typing import Tuple, Dict

import anndata as ad

import config
from embedding.embedding_config import DEFAULT_EMB_COL_NAME


def load_adatas_from_experiment(experiment_name: str, result_dir: Path = config.RESULTS_DIR) -> \
        Tuple[ad.AnnData, Dict[str, ad.AnnData], Dict[str, str]]:
    experiment_dir = Path(result_dir, experiment_name)
    full_adata_file_name = f"adata_with_embedding_in_{DEFAULT_EMB_COL_NAME}.h5ad"
    full_adata_with_emb = ad.read_h5ad(Path(experiment_dir, full_adata_file_name))
    clusters_adatas = {}
    emb_coll_names = {"full": DEFAULT_EMB_COL_NAME}
    for root, subdirectories, _ in os.walk(experiment_dir):
        for sub_dir in [s_dir for s_dir in subdirectories if "cluster" in s_dir]:
            for file in [f for f in os.listdir(Path(root, sub_dir)) if ".h5ad" in f]:
                clusters_adatas[sub_dir] = ad.read_h5ad(Path(root, sub_dir, file))
                emb_coll_names[sub_dir] = file.split(".")[0].split("embedding_in_")[1]
    return full_adata_with_emb, clusters_adatas, emb_coll_names

# example
# full_adata_with_emb, clusters_adatas, emb_coll_names = load_adatas_from_experiment('hierarchical_scvi_2021_11_18__10_10_08')
