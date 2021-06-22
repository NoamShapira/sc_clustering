from pathlib import Path

import anndata as ad
import scanpy as sc
import scvi

import config

result_path = Path(config.RESULTS_DIR, "simple_clustering_2021_06_17__16_16_30", "loaded_data.h5ad")
adata = ad.read(result_path)

sc.pp.filter_genes(adata, min_counts=3)
adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
)

scvi.data.setup_anndata(adata)
model = scvi.model.SCVI(adata)

model.train()

model.save(str(Path(config.RESULTS_DIR, "simple_clustering_2021_06_17__16_16_30")))

