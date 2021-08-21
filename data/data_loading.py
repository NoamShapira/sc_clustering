import logging
from abc import ABC, abstractmethod
from enum import Enum
from functools import partial
from pathlib import Path
from typing import Callable, List, Tuple, Optional

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from tqdm.contrib.concurrent import process_map

import config
from data import meta_data_columns_names
from data.meta_data_columns_names import TREATMENT_ARM
from utils import create_experiment_dir_and_return_path


class AnnDataLoader(ABC):

    @abstractmethod
    def load_data_to_anndata(self) -> ad.AnnData:
        pass

    def load_data_to_anndata_and_save_to_dir(self, results_dir: Optional[Path] = None,
                                             file_name: str = "loaded_data.h5ad") -> Tuple[ad.AnnData, Path]:
        results_dir = results_dir if results_dir is not None else create_experiment_dir_and_return_path(
            "simple_clustering")
        adata = self.load_data_to_anndata()
        adata.write(Path(results_dir, file_name))
        return adata, results_dir


class SeranoDataLoader(AnnDataLoader):
    def __init__(self, experiments_data_dir: Path, meta_data_path: Path,
                 batch_filter_functions: List[Callable[[pd.DataFrame], pd.DataFrame]]):
        self.experiments_data_dir = experiments_data_dir
        self.meta_data_path = meta_data_path
        self.batch_filter_functions = batch_filter_functions if batch_filter_functions is not None else []

    def _get_single_batch(self, row_tpl, col_names):
        row = row_tpl[1]
        cur_data = sc.read_text(Path(self.experiments_data_dir, row[meta_data_columns_names.BATCH_ID] + ".txt"))
        cur_data = cur_data.T
        for col_name in col_names:
            cur_data.obs[col_name] = row[col_name]
        logging.info(f"Reading , batch id - {row[meta_data_columns_names.BATCH_ID]}")
        return cur_data

    def load_data_to_anndata(self) -> ad.AnnData:
        # Read annotation file
        metadata = pd.read_table(self.meta_data_path)
        if config.DEBUG_MODE:
            self.batch_filter_functions.append(lambda df: df.head(config.DEBUG_N_BATCHES))
        for filter_func in self.batch_filter_functions:
            metadata = filter_func(metadata)

        # Read all plates into anndata and merge them
        col_names = metadata.columns
        adatas = process_map(partial(self._get_single_batch, col_names=col_names),
                             list(metadata.iterrows()), max_workers=config.IO_N_WORKERS,
                             desc="loading relevant batches",
                             unit="batch")
        print("merging to single adata")
        adata = ad.concat(adatas, merge="same")
        print(f"converting adata to sparse matrix")
        adata.X = csr_matrix(adata.X)
        print("dropping Mouse columns, some bug with that column")
        adata.obs.drop(['Mouse'], axis='columns', inplace=True)
        return adata


class SeranoDataLoaderDescription(Enum):
    ARMS_1_2_3_FROM_NOAMSH = "arms_1_2_3_from_noamsh"
    ARM_1_FROM_WEINER = "arm_1_from_weiner"


class SeranoDataLoaderFactory:
    @staticmethod
    def _get_only_batches_from_arms_1_2_3(df: pd.DataFrame) -> pd.DataFrame:
        ret_df = df[df[TREATMENT_ARM].aplly(lambda a: a in [1, 2, 3])]
        return ret_df

    @staticmethod
    def create_serano_dataloader(dataloader_description: SeranoDataLoaderDescription) -> SeranoDataLoader:
        if dataloader_description == SeranoDataLoaderDescription.ARMS_1_2_3_FROM_NOAMSH:
            return SeranoDataLoader(experiments_data_dir=config.UMI_DIR_PATH,
                                    meta_data_path=config.UPDATED_META_DATA_PATH,
                                    batch_filter_functions=[SeranoDataLoaderFactory._get_only_batches_from_arms_1_2_3])
        if dataloader_description == SeranoDataLoaderDescription.ARM_1_FROM_WEINER:
            return SeranoDataLoader(experiments_data_dir=config.UMI_DIR_PATH,
                                    meta_data_path=config.META_DATA_PATH,
                                    batch_filter_functions=[])
