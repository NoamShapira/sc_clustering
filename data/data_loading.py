from pathlib import Path
from typing import Callable, List

import anndata as ad
import pandas as pd
import scanpy as sc

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
    adatas = []
    for index, row in metadata.iterrows():
        cur_data = sc.read_text(Path(experiments_data_dir, row[meta_data_columns_names.BATCH_ID] + ".txt"))
        cur_data = cur_data.T
        for col_name in col_names:
            cur_data.obs[col_name] = row[col_name]
        adatas.append(cur_data)
        print(f"Reading , batch id - {row[meta_data_columns_names.BATCH_ID]}, "
              f"batch number - {index}, batch size - {metadata.shape[0]}")

    adata = ad.concat(adatas, merge="same")
    return adata
