import logging
from pathlib import Path

import anndata as ad
import scanpy as sc

import config


def transform_pca_adata(adata: ad.AnnData, pca_svd_solver=config.TL_PCA_SVD_SOLVER, pca_n_comps=config.TL_PCA_N_COMPS):
    sc.tl.pca(adata, svd_solver=pca_svd_solver, n_comps=pca_n_comps)


def compute_neighborhood_graph(adata: ad.AnnData, neighborhood_graph_n_neighbors=config.NEIGHBORHOOD_GRAPH_N_NEIGHBORS,
                               neighborhood_graph_n_pcs=config.NEIGHBORHOOD_GRAPH_N_PCS):
    sc.pp.neighbors(adata, n_neighbors=neighborhood_graph_n_neighbors, n_pcs=neighborhood_graph_n_pcs)
    sc.tl.umap(adata)


def cluster_adata(adata: ad.AnnData, method: str = "leiden", resolution=config.TL_LEIDEN_RESOLUTION):
    adata_copy = adata.copy()
    if method == "leiden":
        sc.tl.leiden(adata, resolution=resolution)
    # if method == "louvain":
    #     sc.tl.louvain(adata, resolution=resolution)
    else:
        raise NotImplementedError
    return adata_copy


def regress_out_and_scale(adata):
    # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
    # Scale the data to unit variance
    sc.pp.regress_out(adata, config.PP_REGRESS_OUT_OBS_KEYS)
    sc.pp.scale(adata, max_value=config.PP_SCALE_MAX_VALUE)


def choose_variable_genes(adata) -> ad.AnnData:
    logging.info("save procced data in adata.raw")
    adata.raw = adata
    # choose only viriable genes
    adata = adata[:, adata.var.highly_variable]
    return adata


def compute_variable_genes(adata, use_triku=False):
    logging.info("compute genes variation")
    sc.pp.highly_variable_genes(adata=adata, **config.PP_HIGHLY_VARIATION_PARAMS)


def normelize_data(adata, normelized_reads_per_cell=config.PP_NORMELIZED_READS_PER_CELL):
    logging.info("normalized reads")
    sc.pp.normalize_total(adata, normelized_reads_per_cell)
    sc.pp.log1p(adata)


def pp_drop_genes_and_cells(adata,
                            min_num_genes_per_cell=config.PP_MIN_NUMBER_OF_GENES_PER_CELL,
                            min_num_cells_per_gene=config.PP_MIN_NUMBER_OF_CELLS_PER_GENE,
                            max_num_genes_per_cell=config.PP_MAX_NUMBER_OF_GENES_PER_CELL,
                            max_num_mt_genes_pre_cell=config.PP_MAX_NUMBER_OF_MT_GENES_PER_CELL):
    logging.info("drop bad cells")
    sc.pp.filter_cells(adata, min_genes=min_num_genes_per_cell)
    sc.pp.filter_cells(adata, max_genes=max_num_genes_per_cell)
    sc.pp.filter_genes(adata, min_cells=min_num_cells_per_gene)
    adata = adata[adata.obs.pct_counts_mt < max_num_mt_genes_pre_cell, :]
    return adata


def pp_rename_vars_add_mt_metrics(adata):
    adata.var_names_make_unique()
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata=adata, **config.PP_QUALITY_CONTROL_PARAMS)


def normalize_and_choose_genes(adata: ad.AnnData):
    new_adata = adata.copy()
    normelize_data(new_adata)
    compute_variable_genes(new_adata)
    new_adata = choose_variable_genes(new_adata)
    regress_out_and_scale(new_adata)
    return new_adata


def create_clusters(adata: ad.AnnData, results_path_dir: Path):
    transform_pca_adata(adata)
    compute_neighborhood_graph(adata)
    cluster_adata(adata)

    adata.write(Path(results_path_dir, "clustered_data.h5ad"))
