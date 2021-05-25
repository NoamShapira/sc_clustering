import logging
from functools import partial
from pathlib import Path
from typing import Callable, List

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from tqdm.contrib.concurrent import thread_map, process_map

import config
from data import meta_data_columns_names


def drop_treated_batches(df: pd.DataFrame) -> pd.DataFrame:
    raise NotImplemented


def get_all_experiments_in_one_anndata(experiments_data_dir: Path = config.UMI_PATH,
                                       meta_data_path: Path = config.META_DATA_PATH) -> ad.AnnData:
    return get_experiments_in_one_anndata(experiments_data_dir, meta_data_path, [])


def get_only_untreated_experiments_in_one_anndata(experiments_data_dir: Path = config.UMI_PATH,
                                                  meta_data_path: Path = config.META_DATA_PATH) -> ad.AnnData:
    return get_experiments_in_one_anndata(experiments_data_dir, meta_data_path, [drop_treated_batches])


def get_experiments_in_one_anndata(experiments_data_dir: Path, meta_data_path: Path,
                                   batch_filter_functions: List[Callable[[pd.DataFrame], pd.DataFrame]]) -> ad.AnnData:
    # Read annotation file
    metadata = pd.read_table(meta_data_path)
    for filter_func in batch_filter_functions:
        metadata = filter_func(metadata)

    # Read all plates into anndata and merge them
    col_names = metadata.columns
    adatas = process_map(partial(get_single_batch, col_names=col_names, experiments_data_dir=experiments_data_dir),
                        list(metadata.iterrows()), max_workers=10)
    # adatas = []
    # for index, row in metadata.iterrows():
    #     cur_data = get_single_batch(col_names, experiments_data_dir, row)
    #     adatas.append(cur_data)
    #     print(f"Reading , batch id - {row[meta_data_columns_names.BATCH_ID]}, "
    #           f"batch number - {index}, batch size - {metadata.shape[0]}")

    adata = ad.concat(adatas, merge="same")
    adata.X = csr_matrix(adata.X)
    return adata


def get_single_batch(row_tpl, col_names, experiments_data_dir):
    row = row_tpl[1]
    cur_data = sc.read_text(Path(experiments_data_dir, row[meta_data_columns_names.BATCH_ID] + ".txt"))
    cur_data = cur_data.T
    for col_name in col_names:
        cur_data.obs[col_name] = row[col_name]
    logging.info(f"Reading , batch id - {row[meta_data_columns_names.BATCH_ID]}")
    return cur_data
