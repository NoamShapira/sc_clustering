import numpy as np
import streamlit as st

from clustering.meta_cell import MetaCellResultsColumnsNames
from meta_cell_vis_streamlit_funcs import full_pipeline_and_load_meta_cell_cached, visualize_meta_cell, \
    visualize_clusters_assosiated_with_meta_cell

adata = full_pipeline_and_load_meta_cell_cached()
st.title("Meta Cells")
meta_cell = st.selectbox("select a metacell", np.unique(adata.obs[MetaCellResultsColumnsNames().meta_cell]))

visualize_meta_cell(adata, meta_cell)

st.title("Scanpy cluster")
visualize_clusters_assosiated_with_meta_cell(adata, meta_cell)
