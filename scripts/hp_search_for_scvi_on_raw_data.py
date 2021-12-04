import argparse
import logging
import os
import sys
from argparse import Namespace
from pathlib import Path

import anndata as ad
import joblib
import optuna
import scanpy as sc
from sklearn.metrics import silhouette_score

sys.path.append('/home/labs/amit/noamsh/repos/sc_clustering')

from clustering.merge_labels_to_adata import merge_labels_to_adata
from embedding.scvi_pipe import train_scvi_on_adata, DEFAULT_EMB_COL_NAME
import config
from data.meta_data_columns_names import TREATMENT_ARM
from utils import get_now_timestemp_as_string

emb_col_name = DEFAULT_EMB_COL_NAME

# input adata, this adata to train scvi on
#   assumed columns : 'Treatment Arm' as batch, 'metacell' - column of labels
asaf_raw_adata_path = Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "cells.h5ad")

# input adata with labels per meta_cell- this script try to optimize embedding to fit labels.
#   assumed columns in labeled_adata.obs: ["cell_type", "broad_cell_type"]
#   and index as metacell
asaf_labeled_adata_path = Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "scanpy_metacells.h5ad")

experiment_name = f"hp_search_scvi_on_raw_gene_counts_{get_now_timestemp_as_string()}"
default_experiment_results_dir_path = Path(config.RESULTS_DIR, experiment_name)

parser = argparse.ArgumentParser(prog='scvi_hp_search_on_raw_data',
                                 description="creating an optuna hyper parameter search")

parser.add_argument("--load_raw_adata_from", type=str, default=str(asaf_raw_adata_path))
parser.add_argument("--load_labeled_adata_from", type=str, default=str(asaf_labeled_adata_path))
parser.add_argument("--results_dir", type=str, default=str(default_experiment_results_dir_path))

parser.add_argument('--n_trials', type=int, default=50)
parser.add_argument('--calc_and_save_umap_of_every_model', type=bool, default=False)
parser.add_argument('--test_mod', type=bool, default=False)

args = parser.parse_args()

experiment_results_dir_path = args.results_dir
new_adata_path = Path(experiment_results_dir_path, "cells.h5ad")
new_labeled_adata_path = Path(experiment_results_dir_path, "scanpy_metacells.h5ad")

os.mkdir(experiment_results_dir_path)
source_adata = ad.read_h5ad(args.load_raw_adata_from)
adata_with_labels = ad.read_h5ad(args.load_labeled_adata_from)
merge_labels_to_adata(source_adata, labels_df=adata_with_labels.obs, labels_col_names=["cell_type", "broad_cell_type"],
                      col_in_adata_to_merge_by="metacell", cols_in_labels_df_to_merge_by="index",
                      cols_to_validate_not_empty=["broad_cell_type"])
source_adata.write(new_adata_path)

sc.settings.figdir = experiment_results_dir_path

adata_with_best_embedding = None
cur_adata_with_embedding = None
cur_model = None
best_model = None


def evaluate_scvi_on_raw_data(trail):
    n_genes = trail.suggest_int("n_genes", 2000, 6000, step=1000)
    epochs = trail.suggest_int("epochs", 50, 100, step=25)
    n_latent = trail.suggest_int("n_latent", 10, 40, step=10)
    n_layers = trail.suggest_int("n_layers", 1, 2)
    n_hidden = trail.suggest_categorical("n_hidden", [32, 64, 128])
    dropout_rate = trail.suggest_float("dropout_rate", 0.1, 0.7)

    if args.test_mod:
        epochs = 5
        n_genes = 500
        n_hidden = 8
        n_layers = 1
        n_latent = 5

    trail_results_path = Path(experiment_results_dir_path,
                              f"hidden_{n_hidden}_layers_{n_layers}_latent_{n_latent}"
                              f"_n_genes_{n_genes}_dropout_{round(dropout_rate, 3)}_epochs_{epochs}")
    os.mkdir(trail_results_path)
    sc.settings.figdir = trail_results_path

    adata = ad.read_h5ad(new_adata_path)
    if args.test_mod:
        adata = adata[0:500, :]

    scvi_args = Namespace(dropout_rate=dropout_rate, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers)
    adata_with_scvi_emb, model = train_scvi_on_adata(adata, args_for_scvi_model=scvi_args,
                                                     return_adata_with_embedding=True, return_model=True,
                                                     emb_col_name=emb_col_name,
                                                     batch_col_name=TREATMENT_ARM, max_epochs=epochs,
                                                     n_variable_genes=n_genes)

    if args.calc_and_save_umap_of_every_model:
        logging.info(f"computing umap for trail {trail_results_path}")
        sc.pp.neighbors(adata, use_rep=emb_col_name)
        sc.tl.umap(adata, min_dist=0.3)
        sc.pl.umap(adata, color=[TREATMENT_ARM, "cell_type"], save="_with_annotation_on_scvi_embedding")

    global cur_adata_with_embedding
    cur_adata_with_embedding = adata_with_scvi_emb.copy()

    global cur_model
    cur_model = model

    return silhouette_score(adata_with_scvi_emb.obsm[emb_col_name],
                            labels=list(adata_with_scvi_emb.obs["broad_cell_type"]))


## not tested for parallelization - some trials can try to edit the best model simultaneously
def update_best_model_callback(study, trial):
    global adata_with_best_embedding
    global best_model
    if study.best_trial == trial:
        adata_with_best_embedding = cur_adata_with_embedding
        best_model = cur_model
        best_model.save(str(Path(experiment_results_dir_path, "best_model_with_batch/")))
        adata_with_best_embedding.write(Path(experiment_results_dir_path, "adata_with_scvi_of_best_model.h5ad"))


study = optuna.create_study(direction="maximize")
study.optimize(evaluate_scvi_on_raw_data,
               n_trials=args.n_trials if not args.test_mod else 1,
               callbacks=[update_best_model_callback])

study_file_path = Path(experiment_results_dir_path, "joblib_study.pkl")
logging.info(f"saving study {study_file_path}")
joblib.dump(study, study_file_path)
logging.info(f"finished saving")

sc.settings.figdir = experiment_results_dir_path
sc.pp.neighbors(adata_with_best_embedding, use_rep=emb_col_name)
sc.tl.umap(adata_with_best_embedding, min_dist=0.3)
sc.pl.umap(adata_with_best_embedding, color=[TREATMENT_ARM, "cell_type"],
           save="_best_with_annotation_on_scvi_embedding")
