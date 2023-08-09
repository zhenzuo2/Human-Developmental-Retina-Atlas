import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_MG_NRPC.h5ad"
)

mv.velocity_embedding_stream(adata_result)
mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

adata_result.obs.loc[:,"root_cells"] = 0
adata_result.obs.loc[(adata_result.obs.Time=="10w")&(adata_result.obs.Region=="Peripheral"),"root_cells"] = 1

mv.latent_time(adata_result)

adata_result.obs["Days"] = adata_result.obs.Days.astype(float)
#adata_result.obsm["X_umap"] = adata[adata_result.obs.index,:].obsm["X_umap"].toarray()
mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="Days",
    alpha=1,
    size=60,
    legend_fontsize=20,
    legend_loc="on data",
    frameon=False,
    legend_fontoutline=True,
    min_mass=3.9,
    smooth=1,
    max_length=100,
    density=2,
    cmap="viridis",
)

fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_mv.velocity_embedding_stream_Days.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="latent_time",
    alpha=1,
    size=60,
    legend_fontsize=20,
    legend_loc="on data",
    frameon=False,
    legend_fontoutline=True,
    min_mass=3.9,
    smooth=1,
    max_length=100,
    density=2,
    title="",
    cmap="viridis",
)

fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_mv.velocity_embedding_stream_latent_time.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="majorclass",
    alpha=1,
    size=60,
    legend_fontsize=20,
    legend_loc="on data",
    frameon=False,
    legend_fontoutline=True,
    min_mass=3.9,
    smooth=1,
    max_length=100,
    density=2,
)

fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_mv.velocity_embedding_stream_majorclass.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)