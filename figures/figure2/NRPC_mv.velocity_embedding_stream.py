import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)

adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/NRPC.h5ad")

mv.velocity_embedding_stream(adata)
adata.obsm['X_umap']=adata_result[adata.obs.index].obsm['X_umap']

adata.obs["Days"] = adata.obs["Days"].astype(float)
mv.velocity_embedding_stream(
    adata,
    basis="umap",
    title = "",
    color="Days",
    alpha=1,
    size=20,
    legend_fontsize=20,
    legend_loc="on data",
    frameon=False,
    legend_fontoutline=True,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/NRPC_velocity_embedding_stream_Days.png",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()