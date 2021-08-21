from pathlib import Path
from typing import NamedTuple

import anndata as ad
import pandas as pd

import config
from clustering.scanpy_cluster import run_full_pipe_from_config
from data.data_loading import SeranoDataLoaderDescription, SeranoDataLoaderFactory


class MetaCellResultsColumnsNames(NamedTuple):
    meta_cell: str = "mc.mc"
    group: str = "group"
    sub_group: str = "sub_group"


class MetaCellResults:
    columns: NamedTuple = MetaCellResultsColumnsNames()

    def __init__(self, df: pd.DataFrame):
        for col_name in self.columns:
            assert col_name in df.columns, f"{col_name} not in dataframe columns"
        self.df: pd.DataFrame = df

    @property
    def obs_names(self):
        return self.df.index


def load_meta_cell_and_merge_to_adata(adata: ad.AnnData, path_to_meta_cell_results: Path) -> ad.AnnData:
    mc_prediction = pd.read_csv(path_to_meta_cell_results)
    mc_prediction = mc_prediction.set_index("Unnamed: 0", drop=True)

    mc_results_prediction = MetaCellResults(mc_prediction)

    ind_of_adata_in_mc = [obs_name in mc_results_prediction.obs_names for obs_name in adata.obs_names]
    combined_adata = adata[ind_of_adata_in_mc, :]

    ind_of_mc_in_adata = [obs_name in adata.obs_names for obs_name in mc_results_prediction.obs_names]
    mc_results_prediction.df = mc_results_prediction.df[ind_of_mc_in_adata]
    for col_name in list(MetaCellResultsColumnsNames()):
        combined_adata.obs[col_name] = mc_results_prediction.df[col_name].fillna("NO_GROUP")
    return combined_adata


def run_full_pipeline_and_load_meta_cell(
        data_loader_desc: SeranoDataLoaderDescription = SeranoDataLoaderDescription.ARM_1_FROM_WEINER,
        path_to_meta_cell_results=config.META_CELL_PATH) -> ad.AnnData:
    raw_adata = SeranoDataLoaderFactory.create_serano_dataloader(data_loader_desc).load_data_to_anndata()
    clustered_adata = run_full_pipe_from_config(raw_adata.copy(), filter_cells_only_during_pp=False)
    return load_meta_cell_and_merge_to_adata(clustered_adata, path_to_meta_cell_results)
