import argparse
import logging
import os
import shutil
import sys
from pathlib import Path

import anndata as ad
import joblib
import optuna
import pandas as pd
import scanpy as sc
import scvi
from sklearn.metrics import silhouette_score

sys.path.append('/home/labs/amit/noamsh/repos/sc_clustering')

import config
from data.meta_data_columns_names import TREATMENT_ARM
from utils import get_now_timestemp_as_string


experiment_name = f"hp_search_scvi_on_raw_gene_counts_{get_now_timestemp_as_string()}"
experiment_results_dir_path = Path(config.RESULTS_DIR, experiment_name)
asaf_raw_adata_path = Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "cells.h5ad")
new_adata_path = Path(experiment_results_dir_path, "cells.h5ad")

parser = argparse.ArgumentParser(prog='scvi_hp_search_on_raw_data',
                                 description="creating an optuna hyper parameter search")
parser.add_argument("--load_raw_adata_from", type=str, default=str(asaf_raw_adata_path))
parser.add_argument('--save_raw_adata_to', type=str, default=str(new_adata_path))
parser.add_argument("--load_labeled_adata_from", type=str,
                    default=str(Path(config.ASAF_META_CELL_ATLAS_DIR_PATH, "scanpy_metacells.h5ad")))
parser.add_argument("--save_labeled_adata_to", type=str,
                    default=str(Path(experiment_results_dir_path, "scanpy_metacells.h5ad")))
parser.add_argument('--n_trials', type=int, default=50)
parser.add_argument('--calc_and_save_umap', type=bool, default=False)
parser.add_argument('--test_mod', type=bool, default=False)

args = parser.parse_args()

os.mkdir(experiment_results_dir_path)
shutil.copyfile(args.load_raw_adata_from, args.save_raw_adata_to)
shutil.copyfile(args.load_labeled_adata_from, args.save_labeled_adata_to)
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

    adata = ad.read_h5ad(args.save_raw_adata_to)
    if args.test_mod:
        adata = adata[0:500, :]
    adata.layers["counts"] = adata.X.copy()  # preserve counts
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_genes,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        batch_key=TREATMENT_ARM
    )

    adata_with_scanpy = ad.read_h5ad(args.save_labeled_adata_to)
    adata.obs = pd.merge(left=adata.obs.astype({"metacell": "str"}).reset_index(),
                         right=adata_with_scanpy.obs[["cell_type", "broad_cell_type"]].reset_index().rename(
                             columns={"index": "mc_num"}),
                         left_on="metacell", right_on="mc_num", how="left", validate="m:1").set_index("index")
    adata_without_nan = adata[~ adata.obs["broad_cell_type"].isna(), :]

    scvi.data.setup_anndata(adata, layer="counts", batch_key=TREATMENT_ARM)
    model = scvi.model.SCVI(adata, n_layers=n_layers, n_latent=n_latent, dropout_rate=dropout_rate, n_hidden=n_hidden)
    model.train(max_epochs=epochs)
    adata.obsm["X_scVI"] = model.get_latent_representation()

    global cur_adata_with_embedding
    cur_adata_with_embedding = adata.copy()

    global cur_model
    cur_model = model

    if args.calc_and_save_umap:
        logging.info(f"computing umap for trail {trail_results_path}")
        sc.pp.neighbors(adata, use_rep="X_scVI")
        sc.tl.umap(adata, min_dist=0.3)
        sc.pl.umap(adata, color=[TREATMENT_ARM, "cell_type"], save="_with_annotation_on_scvi_embedding")

    return silhouette_score(adata_without_nan.obsm["X_scVI"], labels=list(adata_without_nan.obs["broad_cell_type"]))


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
sc.pp.neighbors(adata_with_best_embedding, use_rep="X_scVI")
sc.tl.umap(adata_with_best_embedding, min_dist=0.3)
sc.pl.umap(adata_with_best_embedding, color=[TREATMENT_ARM, "cell_type"],
           save="_best_with_annotation_on_scvi_embedding")
