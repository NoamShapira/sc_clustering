from typing import List

import anndata as ad
import pandas as pd


def merge_labels_to_adata(adata: ad.AnnData, labels_df: pd.DataFrame, labels_col_names: List[str],
                           col_in_adata_to_merge_by: str, cols_in_labels_df_to_merge_by: str,
                           cols_to_validate_not_empty: List[str]):
    if cols_in_labels_df_to_merge_by == "index":
        labels_df = labels_df.reset_index()
    adata.obs = pd.merge(
        left=adata.obs.astype({col_in_adata_to_merge_by: "str"}).reset_index(),
        right=labels_df[labels_col_names], how="left", validate="m:1",
        left_on=col_in_adata_to_merge_by, right_on=cols_in_labels_df_to_merge_by).set_index("index")
    adata = adata[~ adata.obs[cols_to_validate_not_empty].isna().any(axis=1), :]