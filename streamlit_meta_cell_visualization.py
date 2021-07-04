import streamlit as st
import numpy as np

from clustering.meta_cell import meta_cell_and_clustering_comparison, MetaCellResultsColumnsNames

adata = meta_cell_and_clustering_comparison()
meta_cell = st.selectbox("select a metacell", np.unipue(adata.obs[MetaCellResultsColumnsNames().meta_cell]))
visualize_meta_cell(adata, meta_cell)
visualize_clusters_assosiated_with_meta_cell(adata, meta_cell)

