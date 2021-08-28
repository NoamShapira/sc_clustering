import pandas as pd

import config
from data.data_loaders import SeranoDataLoader
from data.meta_data_columns_names import TREATMENT_ARM
from data.serono_data_loading_description import SeranoDataLoaderDescription


class SeranoDataLoaderFactory:
    @staticmethod
    def _get_only_batches_from_arms_1_2_3(df: pd.DataFrame) -> pd.DataFrame:
        ret_df = df[df[TREATMENT_ARM].apply(lambda a: a in [1, 2, 3])]
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
