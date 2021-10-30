from argparse import Namespace
from pathlib import Path
from typing import Optional

import anndata as ad
import scanpy as sc
import scvi


def train_scvi_on_adata(adata: ad.AnnData, n_variable_genes: int, batch_col_name: str,
                        return_adata_with_embedding: bool, emb_col_name: str, args_for_scvi_model: Namespace,
                        norm_target_sum: float = 1e6, results_dir: Optional[Path] = None,
                        raw_counts_to_layer_name: str = "counts",
                        max_epochs: Optional[int] = None, drop_genes_and_cells: bool = False) -> Optional[ad.AnnData]:
    if drop_genes_and_cells:
        raise NotImplementedError

    adata.layers[raw_counts_to_layer_name] = adata.X.copy()  # preserve counts
    sc.pp.normalize_total(adata, target_sum=norm_target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_variable_genes,
        subset=True,
        layer=raw_counts_to_layer_name,
        flavor="seurat_v3",
        batch_key=batch_col_name
    )

    scvi.data.setup_anndata(adata, layer=raw_counts_to_layer_name, batch_key=batch_col_name)
    model = scvi.model.SCVI(adata, n_layers=args_for_scvi_model.n_layers, n_latent=args_for_scvi_model.n_latent,
                            dropout_rate=args_for_scvi_model.dropout_rate, n_hidden=args_for_scvi_model.n_hidden)
    model.train(max_epochs=max_epochs)
    adata.obsm[emb_col_name] = model.get_latent_representation()

    if results_dir is not None:
        adata.write(Path(results_dir, f"adata_with_embedding_in_{emb_col_name}.h5ad"))
        model.save(str(Path(results_dir, "model/")))

    if return_adata_with_embedding:
        return adata.copy()
