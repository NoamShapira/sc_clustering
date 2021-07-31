import logging
from pathlib import Path

import anndata as ad
import scanpy as sc

import config
from data.data_loading import load_data_and_save_to_results_dir
from utils import get_now_timestemp_as_string


def transform_pca_adata(adata: ad.AnnData, pca_svd_solver=config.TL_PCA_SVD_SOLVER, pca_n_comps=config.TL_PCA_N_COMPS):
    sc.tl.pca(adata, svd_solver=pca_svd_solver, n_comps=pca_n_comps)


def compute_neighborhood_graph(adata: ad.AnnData, neighborhood_graph_n_neighbors=config.NEIGHBORHOOD_GRAPH_N_NEIGHBORS,
                               neighborhood_graph_n_pcs=config.NEIGHBORHOOD_GRAPH_N_PCS):
    sc.pp.neighbors(adata, n_neighbors=neighborhood_graph_n_neighbors, n_pcs=neighborhood_graph_n_pcs)
    sc.tl.umap(adata)


def cluster_adata(adata: ad.AnnData, method: str = "leiden", resolution=config.TL_LEIDEN_RESOLUTION):
    if method == "leiden":
        sc.tl.leiden(adata, resolution=resolution)
    # if method == "louvain":
    #     sc.tl.louvain(adata, resolution=resolution)
    else:
        raise NotImplementedError


def choose_variable_genes(adata) -> ad.AnnData:
    logging.info("save procced data in adata.raw")
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    return adata


def compute_variable_genes(adata, use_triku=False):
    logging.info("compute genes variation")
    sc.pp.highly_variable_genes(adata=adata, **config.PP_HIGHLY_VARIATION_PARAMS)


def normelize_data(adata, normelized_reads_per_cell=config.PP_NORMELIZED_READS_PER_CELL):
    logging.info("normalized reads")
    sc.pp.normalize_total(adata, normelized_reads_per_cell)
    sc.pp.log1p(adata)


def pp_drop_genes_and_cells(adata, filter_cells_only_during_pp=False,
                            min_num_genes_per_cell=config.PP_MIN_NUMBER_OF_GENES_PER_CELL,
                            min_num_cells_per_gene=config.PP_MIN_NUMBER_OF_CELLS_PER_GENE,
                            max_num_genes_per_cell=config.PP_MAX_NUMBER_OF_GENES_PER_CELL,
                            max_num_mt_genes_pre_cell=config.PP_MAX_PCT_OF_MT_GENES_PER_CELL,
                            min_num_counts_per_cell=config.PP_MIN_NUMBER_OF_GENES_PER_CELL):
    logging.info("drop bad cells")
    sc.pp.filter_cells(adata, min_genes=min_num_genes_per_cell)
    sc.pp.filter_cells(adata, max_genes=max_num_genes_per_cell)
    sc.pp.filter_cells(adata, min_counts=min_num_counts_per_cell)
    adata = adata[adata.obs.pct_counts_mt < max_num_mt_genes_pre_cell, :]
    if not filter_cells_only_during_pp:
        logging.info("drop bad genes")
        sc.pp.filter_genes(adata, min_cells=min_num_cells_per_gene)
    return adata


def pp_rename_vars_add_mt_metrics(adata):
    adata.var_names_make_unique()
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata=adata, **config.PP_QUALITY_CONTROL_PARAMS)


def normalize_and_choose_genes(adata: ad.AnnData, drop_unvariable_genes: bool = config.PP_DROP_UNVARIABLE_GENES,
                               regress_out_total_cont_and_mt: bool = config.PP_REGRESS_OUT_TOTAL_CONT_AND_MT):
    normelize_data(adata)
    if drop_unvariable_genes:
        compute_variable_genes(adata)
        adata = choose_variable_genes(adata)
    if regress_out_total_cont_and_mt:
        sc.pp.regress_out(adata, config.PP_REGRESS_OUT_OBS_KEYS)
    sc.pp.scale(adata, max_value=config.PP_SCALE_MAX_VALUE)


def run_full_pipe_from_config(adata: ad.AnnData, filter_cells_only_during_pp: bool):
    adata = pp_choose_genes_and_normelize(adata, filter_cells_only_during_pp)
    transform_pca_adata(adata)
    compute_neighborhood_graph(adata)
    cluster_adata(adata)
    return adata


def pp_choose_genes_and_normelize(adata, filter_cells_only_during_pp: bool = False):
    pp_rename_vars_add_mt_metrics(adata)
    adata = pp_drop_genes_and_cells(adata, filter_cells_only_during_pp)
    normalize_and_choose_genes(adata)
    return adata


def run_full_clustering_pipe_and_create_results_dir():
    adata, experiment_results_dir_path = load_data_and_save_to_results_dir()
    adata = run_full_pipe_from_config(adata, filter_cells_only_during_pp=True)
    adata.write(Path(experiment_results_dir_path, f"final_adata_{get_now_timestemp_as_string()}.h5ad"))
