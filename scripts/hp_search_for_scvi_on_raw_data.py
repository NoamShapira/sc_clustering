import argparse
import os
import shutil
from pathlib import Path

import anndata as ad
import joblib
import optuna
import scanpy as sc
import pandas as pd
import scvi
from sklearn.metrics import silhouette_score

import config
from data.meta_data_columns_names import TREATMENT_ARM
from utils import get_now_timestemp_as_string

experiment_name = f"scvi_on_raw_gene_counts_{get_now_timestemp_as_string()}"
experiment_results_dir_path = Path(config.RESULTS_DIR, experiment_name)
os.mkdir(experiment_results_dir_path)

asaf_meta_cell_dir_path = Path("/home/labs/amit/weiner/Serono/Serono14/clustering_results/simple_clustering_2021_09_05__15_04_43")
old_adata_path = Path(asaf_meta_cell_dir_path, "cells.h5ad")
new_adata_path = Path(experiment_results_dir_path, "cells.h5ad")

parser = argparse.ArgumentParser(prog='scvi_hp_search_on_raw_data',
                                 description="creating an optuna hyper parameter search")
parser.add_argument("--raw_adata_path", type=str, default=str(old_adata_path))
parser.add_argument('--new_adata_path', type=str, default=str(new_adata_path))
parser.add_argument('--n_trials', type=int, default=50)

args = parser.parse_args()

shutil.copyfile(args.old_adata_path, args.new_adata_path)
sc.settings.figdir = experiment_results_dir_path

def evaluate_scvi_on_raw_data(trail):
    n_genes = 2000
    epochs = trail.suggest_int("epochs", 50, 200, step=50)
    n_latent  = trail.suggest_int("n_latent", 10, 50, step=10)
    n_layers = trail.suggest_int("n_layers", 1, 2)
    n_hidden = trail.suggest_categorical("n_hidden", [32, 64, 128, 256])
    dropout_rate = trail.suggest_float("dropout_rate", 0.1, 0.7)

    adata = ad.read_h5ad(new_adata_path)
    adata.layers["counts"] = adata.X.copy() # preserve counts
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
    scvi.data.setup_anndata(adata, layer="counts", batch_key=TREATMENT_ARM)

    model = scvi.model.SCVI(adata, n_layers=n_layers, n_latent=n_latent, dropout_rate=dropout_rate, n_hidden=n_hidden)

    model.train(epochs=epochs)
    model.save(str(Path(experiment_results_dir_path, "model_with_batch/")))
    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.write(Path(experiment_results_dir_path, "adata_with_scvi.h5ad"))

    adata_with_scanpy = ad.read_h5ad(Path(asaf_meta_cell_dir_path, "scanpy_metacells.h5ad"))
    adata.obs = pd.merge(left=adata.obs,
                         right=adata_with_scanpy.obs["cell_type"].reset_index().astype({"index": 'int32'}),
                         left_on="metacell", right_on="index",
                         how="left")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata, min_dist=0.3)
    sc.pl.umap(adata, color=[TREATMENT_ARM, "cell_type"])

    adata_without_nan = adata[~ adata.obs["broad_cell_type"].isna(), :]
    return silhouette_score(adata_without_nan.obsm["X_scVI"], labels=list(adata_without_nan.obs["broad_cell_type"]))

study = optuna.create_study(direction="maximize") # # sampler=optuna.samplers.GridSampler(search_space)
study.optimize(evaluate_scvi_on_raw_data, n_trials=args.n_trials)

study_file_path = Path(experiment_results_dir_path, "joblib_study.pkl")
print(f"saving study {study_file_path}")
joblib.dump(study, study_file_path)

