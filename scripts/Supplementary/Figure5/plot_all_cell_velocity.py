import os
import scvelo as scv
import joblib
import cellrank as cr
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad",
    cache=False,
)

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=50000)
scv.pp.log1p(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=50000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    legend_fontsize=20,
    title="",
    color="majorclass",
    size=10,
    colorbar=True,
    dpi=600,
    figsize=(10, 10),
    linewidth=5,
    alpha=1,
    legend_loc="right margin",
    #arrow_color = "white",
    palette={
    "MG": "#9467bd",
    "Rod": "#17becf",
    "BC": "#bcbd22",
    "RGC": "#d62728",
    "NRPC": "#ff7f0e",
    "Cone": "#e377c2",
    "HC": "#2ca02c",
    "PRPC": "#1f77b4",
    "AC": "#8c564b",
}
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure5/scv.pl.velocity_embedding_stream_scpred_prediction.tiff",
    bbox_inches="tight",
    transparent=True,
)