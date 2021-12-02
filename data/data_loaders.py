import logging
import os
from abc import ABC, abstractmethod
from functools import partial
from pathlib import Path
from typing import Callable, List, Optional

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from tqdm.contrib.concurrent import process_map

from data.meta_data_columns_names import BATCH_ID


class AnnDataLoader(ABC):

    @abstractmethod
    def load_data_to_anndata(self) -> ad.AnnData:
        pass


DEBUG_MODE = False
DEBUG_N_BATCHES = 10

class AmpBatchDataLoader(AnnDataLoader):
    def __init__(self, experiments_data_dir: Path, meta_data_path: Path,
                 batch_filter_functions: Optional[List[Callable[[pd.DataFrame], pd.DataFrame]]] = None):
        self.experiments_data_dir = experiments_data_dir
        self.meta_data_path = meta_data_path
        self.batch_filter_functions = batch_filter_functions if batch_filter_functions is not None else []

    def _get_single_batch(self, row_tpl, col_names):
        row = row_tpl[1]
        cur_data = sc.read_text(Path(self.experiments_data_dir, row[BATCH_ID] + ".txt"))
        cur_data = cur_data.T
        for col_name in col_names:
            cur_data.obs[col_name] = row[col_name]
        logging.info(f"Reading , batch id - {row[BATCH_ID]}")
        cur_data.X = csr_matrix(cur_data.X)
        logging.info(f"converting bath batch id - {row[BATCH_ID]} to sparse matrix")
        return cur_data

    def load_data_to_anndata(self) -> ad.AnnData:
        # Read annotation file
        if self.meta_data_path.suffix == "csv":
            metadata_df = pd.read_table(self.meta_data_path)
        else:
            metadata_df = pd.read_excel(self.meta_data_path)
        metadata_df = metadata_df.dropna(axis=0, subset=[BATCH_ID])
        if DEBUG_MODE:
            self.batch_filter_functions.append(lambda df: df.head(DEBUG_N_BATCHES))
        for filter_func in self.batch_filter_functions:
            metadata_df = filter_func(metadata_df)

        # Read all plates into anndata and merge them
        col_names = metadata_df.columns
        if DEBUG_MODE:
            adatas = [self._get_single_batch(batch_tpl, col_names) for batch_tpl in metadata_df.iterrows()]
        else:
            adatas = process_map(partial(self._get_single_batch, col_names=col_names),
                                 list(metadata_df.iterrows()), max_workers=os.cpu_count()//2,
                                 desc="loading relevant batches",
                                 unit="batch")
        logging.info("merging to single adata")
        adata = ad.concat(adatas, merge="same")

        bad_columns = ['Mouse'] #, "Experimental Batch", 'Libraries date','Frozen']
        adata.obs = adata.obs.astype(str)
        logging.info(f"dropping - {str(bad_columns)} columns, some bug with that columns")
        adata.obs.drop(bad_columns, axis='columns', inplace=True)
        return adata
