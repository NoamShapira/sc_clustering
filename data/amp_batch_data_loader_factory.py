from enum import Enum

import pandas as pd

import config
from data.data_loaders import AmpBatchDataLoader
from data.meta_data_columns_names import TREATMENT_ARM


class PlatesLoaderDescription(Enum):
    ARMS_1_2_3_FROM_NOAMSH = "arms_1_2_3_from_noamsh"
    ARM_1_FROM_WEINER = "arm_1_from_weiner"


class PlatesDataLoaderFactory:
    @staticmethod
    def _get_only_batches_from_arms_1_2_3(df: pd.DataFrame) -> pd.DataFrame:
        ret_df = df[df[TREATMENT_ARM].apply(lambda a: a in [1, 2, 3])]
        return ret_df

    @staticmethod
    def create_amp_batch_dataloader(dataloader_description: PlatesLoaderDescription) -> AmpBatchDataLoader:
        if dataloader_description == PlatesLoaderDescription.ARMS_1_2_3_FROM_NOAMSH:
            return AmpBatchDataLoader(experiments_data_dir=config.UMI_DIR_PATH,
                                      meta_data_path=config.UPDATED_META_DATA_PATH,
                                      batch_filter_functions=[PlatesDataLoaderFactory._get_only_batches_from_arms_1_2_3])
        if dataloader_description == PlatesLoaderDescription.ARM_1_FROM_WEINER:
            return AmpBatchDataLoader(experiments_data_dir=config.UMI_DIR_PATH,
                                      meta_data_path=config.META_DATA_PATH,
                                      batch_filter_functions=[])
