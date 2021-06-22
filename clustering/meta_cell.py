from pathlib import Path
from typing import List

import anndata as ad
import pandas as pd


def load_meta_cell_and_merge_to_adata(adata: ad.AnnData, path_to_meta_cel_Results: Path,
                                      reference_col_names: List[str]) -> ad.AnnData:
    mc_prediction = pd.read_csv(path_to_meta_cel_Results)
    mc_prediction = mc_prediction.set_index("Unnamed: 0", drop=True)

    ind_of_adata_in_mc = [obs_name in mc_prediction.index for obs_name in adata.obs_names]
    combined_adata = adata[ind_of_adata_in_mc, :]

    ind_of_mc_in_adata = [obs_name in adata.obs_names for obs_name in mc_prediction.index]
    mc_prediction = mc_prediction[ind_of_mc_in_adata]
    for col_name in reference_col_names:
        combined_adata.obs[col_name] = mc_prediction[col_name].fillna("NO_GROUP")
    return combined_adata
