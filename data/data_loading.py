import logging
from functools import partial
from pathlib import Path
from typing import Callable, List, Tuple

import anndata
import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from tqdm.contrib.concurrent import process_map

import config
from data import meta_data_columns_names
from data.meta_data_columns_names import TREATMENT_ARM
from utils import create_experiment_dir_and_return_path


def get_only_batches_from_arms_1_2_3(df: pd.DataFrame) -> pd.DataFrame:
    ret_df = df[df[TREATMENT_ARM].aplly(lambda a: a in [1, 2, 3])]
    return ret_df


def get_all_experiments_in_one_anndata(experiments_data_dir: Path = config.UMI_DIR_PATH,
                                       meta_data_path: Path = config.META_DATA_PATH) -> ad.AnnData:
    return get_experiments_in_one_anndata(experiments_data_dir, meta_data_path, [])


def get_arms_1_2_3_in_one_anndata(experiments_data_dir: Path = config.UMI_DIR_PATH,
                                  meta_data_path: Path = config.UPDATED_MEATADATA_PATH) -> ad.AnnData:
    return get_experiments_in_one_anndata(experiments_data_dir, meta_data_path, [get_only_batches_from_arms_1_2_3])


def get_experiments_in_one_anndata(experiments_data_dir: Path, meta_data_path: Path,
                                   batch_filter_functions: List[Callable[[pd.DataFrame], pd.DataFrame]]) -> ad.AnnData:
    # Read annotation file
    metadata = pd.read_table(meta_data_path)
    if config.DEBUG_MODE:
        batch_filter_functions.append(lambda df: df.head(config.DEBUG_N_BATCHES))
    for filter_func in batch_filter_functions:
        metadata = filter_func(metadata)

    # Read all plates into anndata and merge them
    col_names = metadata.columns
    adatas = process_map(partial(get_single_batch, col_names=col_names, experiments_data_dir=experiments_data_dir),
                         list(metadata.iterrows()), max_workers=config.IO_N_WORKERS, desc="loading relevant batches",
                         unit="batch")
    print("merging to single adata")
    adata = ad.concat(adatas, merge="same")
    print(f"converting adata to sparse matrix")
    adata.X = csr_matrix(adata.X)
    print("dropping Mouse columns, some bug with that column")
    adata.obs.drop(['Mouse'], axis='columns', inplace=True)
    return adata


def get_single_batch(row_tpl, col_names, experiments_data_dir):
    row = row_tpl[1]
    cur_data = sc.read_text(Path(experiments_data_dir, row[meta_data_columns_names.BATCH_ID] + ".txt"))
    cur_data = cur_data.T
    for col_name in col_names:
        cur_data.obs[col_name] = row[col_name]
    logging.info(f"Reading , batch id - {row[meta_data_columns_names.BATCH_ID]}")
    return cur_data


LOADING_FUNCTIONS = {
    "arms_1_2_3_from_noamsh": get_arms_1_2_3_in_one_anndata,
    "arm_1_from_weiner": get_all_experiments_in_one_anndata
}


def load_data_and_save_to_results_dir(loading_func_name: str = None) -> Tuple[anndata.AnnData, Path]:
    experiment_results_dir_path = create_experiment_dir_and_return_path("simple_clustering")
    if loading_func_name is not None:
        assert loading_func_name in LOADING_FUNCTIONS.keys()
    else:
        loading_func_name = "arm_1_from_weiner"
    adata = LOADING_FUNCTIONS[loading_func_name]()
    adata.write(Path(experiment_results_dir_path, "loaded_data.h5ad"))
    return adata, experiment_results_dir_path
