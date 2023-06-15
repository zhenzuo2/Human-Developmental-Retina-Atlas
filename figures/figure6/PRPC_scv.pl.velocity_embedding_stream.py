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
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC.h5ad"
)

mv.velocity_embedding_stream(adata_result)
mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

adata_result.obs["Days"] = adata_result.obs.Days.astype(float)

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
    min_mass=3.9,
    smooth=1,
    max_length=100,
    density=2,
)

fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "PRPC_mv.velocity_embedding_stream_Days.svg",
    bbox_inches="tight",
    backend="cairo",
    transparent=True,
)

mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="latent_time",
    alpha=1,
    size=20,
    legend_fontsize=20,
    legend_loc="on data",
    frameon=False,
    legend_fontoutline=True,
    min_mass=3.9,
    smooth=1,
    max_length=100,
    density=2,
    title = "",cmap = 'coolwarm'
)

fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "PRPC_mv.velocity_embedding_stream_latent_time.svg",
    bbox_inches="tight",
    backend="cairo",
    transparent=True,
)


mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="majorclass",
    alpha=1,
    size=20,
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
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "PRPC_mv.velocity_embedding_stream_majorclass.svg",
    bbox_inches="tight",
    backend="cairo",
    transparent=True,
)