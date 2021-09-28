import os
import sys
from pathlib import Path

import anndata as ad
import joblib
import numpy as np
import optuna
import pytorch_lightning.loggers as pl_loggers
import scanpy as sc
import scvi
from sklearn.metrics import silhouette_score

sys.path.append('/home/labs/amit/noamsh/repos/sc_clustering')

import config
from utils import get_now_timestemp_as_string

metacells_file_path = Path(config.RESULTS_DIR, "mc_2021_09_23__09_59_00", "metacells.h5ad")
metacells_with_annot_file_path = Path(config.RESULTS_DIR, "mc_2021_09_23__09_59_00", "scanpy_metacells.h5ad")

experiment_name = f"hp_search_scvi_on_meta_cells_{get_now_timestemp_as_string()}"
experiment_results_dir_path = Path(config.RESULTS_DIR, experiment_name)
os.mkdir(experiment_results_dir_path)

tb_logs_path = Path(experiment_results_dir_path, "tb_logs")
os.mkdir(tb_logs_path)

normelize_data = False
choose_features_from_mc = True
N_TRIALS = 50

search_space = {
    "n_hidden": [128, 64],
    "n_layers": [1, 2],
    "n_latent": [10, 25],
    "n_neighbors": [10, 15],
    'dropout_rate': [0.1, 0.3],
    'epochs': [400]
}

def objective(trail):
    n_hidden = trail.suggest_categorical("n_hidden", [32, 64, 128, 256])
    n_layers = trail.suggest_int("n_layers", 1, 2)
    n_latent = trail.suggest_int("n_latent", 10, 30, step=5)
    n_neighbors = trail.suggest_int("n_neighbors", 10, 30)
    dropout_rate = trail.suggest_float("dropout_rate", 0.1, 0.7)
    epochs = trail.suggest_int("epochs", 200, 800, step=200)

    trail_results_path = Path(experiment_results_dir_path,
                              f"hidden_{n_hidden}_layers_{n_layers}_latent_{n_latent}"
                              f"_neighbors_{n_neighbors}_dropout_{dropout_rate}_epochs_{epochs}")
    os.mkdir(trail_results_path)
    sc.settings.figdir = trail_results_path

    adata = ad.read_h5ad(metacells_file_path)
    adata_with_scanpy = ad.read_h5ad(metacells_with_annot_file_path)
    adata.obs['broad_cell_type'] = adata_with_scanpy.obs['broad_cell_type']
    adata.obs['leiden'] = adata_with_scanpy.obs['leiden']
    adata.obs['cell_type'] = adata_with_scanpy.obs['cell_type']

    sc.pp.filter_cells(adata, min_genes=200)  # filter cells with fewer than 200 genes
    sc.pp.filter_cells(adata,
                       min_counts=200)  # This is a weaker threshold than above. It is just to population the n_counts column in adata
    sc.pp.filter_genes(adata, min_cells=2)  # filter genes detected in fewer than 3 cells

    if normelize_data:
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.1, max_mean=3, min_disp=0.8)
        adata.var['highly_variable'] = np.logical_or(adata.var['highly_variable'], adata.var.feature_gene > 3)
    if choose_features_from_mc:
        adata = adata[:,
                [gn for gn in adata.var_names if
                 gn in adata_with_scanpy.var['feature_gene'] and adata_with_scanpy.var['feature_gene'][gn] != 0]].copy()

    scvi.data.setup_anndata(adata)
    model = scvi.model.SCVI(adata, n_hidden=n_hidden, dropout_rate=dropout_rate, n_layers=n_layers, n_latent=n_latent)

    tb_logger = pl_loggers.TensorBoardLogger(str(Path(trail_results_path, "logs/")))
    model.train(max_epochs=epochs, logger=tb_logger)
    adata.obsm["X_scVI"] = model.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=n_neighbors)
    sc.tl.umap(adata, min_dist=0.2)
    sc.pl.umap(adata, color=["cell_type", "leiden", "broad_cell_type"],
               save="_annotation_on_graph_with_scvi_embedding")

    adata.write(Path(trail_results_path, "adata_with_annot_and_scvi.h5ad"))

    return silhouette_score(adata.obsm["X_scVI"], labels=list(adata.obs["broad_cell_type"]))


study = optuna.create_study(direction="maximize") # # sampler=optuna.samplers.GridSampler(search_space)
study.optimize(objective, n_trials=N_TRIALS)

joblib.dump(study, Path(experiment_results_dir_path, "study.pkl"))