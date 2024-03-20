import scanpy as sc
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import os
import scipy
import numpy as np
import scvelo as scv
import multivelo as mv

mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["scatter.edgecolors"] = "none"
font = {"family": "Arial", "size": 12}
plt.rc("font", **font)

adata_result = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/NRPC.h5ad"
)
adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/NRPC_res2.h5ad"
)

adata_result.obsm["X_umap"] = adata[adata_result.obs.index].obsm["X_umap"]
adata_result.obs["subclass"] = adata[adata_result.obs.index].obs["subclass"]

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)

mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="Days",
    alpha=1,
    size=20,
    legend_fontsize=20,
    legend_loc="on data",
    frameon=False,
    legend_fontoutline=True,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_Days.tiff",
    transparent=True,
)
