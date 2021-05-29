import logging
from pathlib import Path

import anndata as ad
import scanpy as sc

import config


def pre_procces_adata(adata: ad.AnnData):
    adata = pp_drop_genes_and_cells(adata)

    normelize_data(adata)

    compute_variable_genes(adata)

    adata = choose_variable_genes(adata)

    # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
    # Scale the data to unit variance
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)


def choose_variable_genes(adata):
    logging.info("save procced data in adata.raw")
    adata.raw = adata
    # choose only viriable genes
    adata = adata[:, adata.var.highly_variable]
    return adata


def compute_variable_genes(adata):
    logging.info("compute genes variation")
    sc.pp.highly_variable_genes(adata=adata, **config.PP_HIGHLY_VARIATION_PARAMS)


def normelize_data(adata):
    logging.info("normalized reads")
    sc.pp.normalize_total(adata, config.PP_NORMELIZED_READS_PER_CELL)
    sc.pp.log1p(adata)


def pp_drop_genes_and_cells(adata):
    logging.info("drop bad cells")
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=config.PP_MIN_NUMBER_OF_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=config.PP_MIN_NUMBER_OF_CELLS_PER_GENE)
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata=adata, **config.PP_QUALITY_CONTROL_PARAMS)
    adata = adata[adata.obs.n_genes_by_counts < config.PP_MAX_NUMBER_OF_GENES_PER_CELL, :]
    adata = adata[adata.obs.pct_counts_mt < config.PP_MAX_NUMBER_OF_MT_GENES_PER_CELL, :]
    return adata


def create_clusters(adata: ad.AnnData, results_path_dir: Path):
    raise NotImplementedError
