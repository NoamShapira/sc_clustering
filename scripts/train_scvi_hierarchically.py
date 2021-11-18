import argparse
import logging
import os
import shutil
import sys
from pathlib import Path

import anndata as ad

sys.path.append('/home/labs/amit/noamsh/repos/sc_clustering')

from clustering.split_adata_to_major_clusters import split_to_major_cluster
import config
from data.meta_data_columns_names import TREATMENT_ARM
from embedding.scvi_pipe import train_scvi_on_adata
from utils import get_now_timestemp_as_string

asaf_raw_adata_path = Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "cells.h5ad")
asaf_labeled_adata_path = Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "scanpy_metacells.h5ad")

experiment_name = f"hierarchical_scvi_{get_now_timestemp_as_string()}"
experiment_results_dir_path = Path(config.RESULTS_DIR, experiment_name)

parser = argparse.ArgumentParser()
parser.add_argument('--test_mod', type=bool, default=False)

parser.add_argument("--load_raw_adata_from", type=str, default=str(asaf_raw_adata_path))
parser.add_argument("--load_labeled_adata_from", type=str, default=str(asaf_labeled_adata_path))
parser.add_argument("--save_results_to", type=str, default=str(experiment_results_dir_path))

# spliting params
parser.add_argument("--num_clusters_to_split_to", type=int, default=4)
parser.add_argument("--clustering_method", type=str, default="leiden_and_agglomerative")
parser.add_argument("--min_percentage_of_cluster_to_compute_new_embedding", type=float, default=0.1)
parser.add_argument("--split_and_continue_from_dir", type=str, default="")

# preprocessing
parser.add_argument("--pp_n_top_variable_genes", type=int, default=3000)

# scvi params
parser.add_argument("--n_layers", type=int, default=1)
parser.add_argument("--n_latent", type=int, default=30)
parser.add_argument("--dropout_rate", type=float, default=0.5)
parser.add_argument("--n_hidden", type=int, default=64)

args = parser.parse_args()

source_adata_path = Path(args.save_results_to, "cells.h5ad")
source_labeled_adata_path = Path(args.save_results_to, "scanpy_metacells.h5ad")

if os.path.isdir(args.save_results_to):
    raise ValueError(f"an existing results dir was passed as arg,"
                     f" expect non existing dir and the script will create it")
os.mkdir(args.save_results_to)
shutil.copyfile(args.load_raw_adata_from, source_adata_path)
shutil.copyfile(args.load_labeled_adata_from, source_labeled_adata_path)

adata = ad.read_h5ad(source_adata_path)
# adata_with_labels = ad.read_h5ad(source_labeled_adata_path)
if args.test_mod:
    adata = adata[0:300, :]
    # adata_with_labels = adata_with_labels[adata.obs, :]

emb_col_name = "X_scVI"

if args.split_and_continue_from_dir == "":
    adata_with_scvi_emb = train_scvi_on_adata(adata, args_for_scvi_model=args, results_dir=args.save_results_to,
                                              return_adata_with_embedding=True, emb_col_name=emb_col_name,
                                              batch_col_name=TREATMENT_ARM,
                                              n_variable_genes=args.pp_n_top_variable_genes)
else:
    logging.info("training scvi model on all data")
    adata_with_scvi_emb = ad.read_h5ad(
        Path(args.split_and_continue_from_dir, f"adata_with_embedding_in_{emb_col_name}.h5ad"))
    adata_with_scvi_emb.write(Path(args.save_results_to, f"loaded_adata_with_embedding_in_{emb_col_name}.h5ad"))

logging.info("spliting adata to major clusters")
adatas_of_major_cluster = split_to_major_cluster(adata_with_scvi_emb, emb_col_name=emb_col_name,
                                                 clustering_model_name=args.clustering_method,
                                                 num_clusters=args.num_clusters_to_split_to)

for i, adata_of_cluster in enumerate(adatas_of_major_cluster):
    split_results_dir_path = Path(args.save_results_to, f"cluster_{i}")
    os.mkdir(split_results_dir_path)
    if len(adata_of_cluster) >= args.min_percentage_of_cluster_to_compute_new_embedding * len(adata_with_scvi_emb):
        logging.info(f"training scvi model on cluster {i}")
        train_scvi_on_adata(adata_of_cluster, args_for_scvi_model=args,
                            results_dir=split_results_dir_path, return_adata_with_embedding=False,
                            emb_col_name=f"{emb_col_name}_{i}", batch_col_name=TREATMENT_ARM,
                            n_variable_genes=args.pp_n_top_variable_genes)
    else:
        adata_of_cluster.write(Path(split_results_dir_path, f"adata_with_embedding_in_{emb_col_name}.h5ad"))
        logging_msg = f"for cluster {i}, did not compute embedding. the cluster was too small"
        logging.info(logging_msg)
        print(logging_msg)

print("finished script !!")
