import scvelo as scv
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/"

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata_dynamics.h5ad"
)
adata.obs["Days"] = adata.obs["Days"].astype(float)

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)
adata.obs["scpred_prediction"] = pd.Categorical(
    list(adata.obs["scpred_prediction"]),
    categories=[
        "RPC",
        "RGC",
        "Cone",
        "HC",
        "AC",
        "Rod",
        "BC",
        "MG",
    ],
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    legend_fontsize=20,
    title="",
    color="Days",
    size=10,
    colorbar=True,
    dpi=600,
    figsize=(10, 10),
    linewidth=5,
    alpha=1,
    legend_loc="right margin",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "woadult_scv.pl.velocity_embedding_stream_Days.svg",
    bbox_inches="tight",
    backend = "cairo",
    transparent=True,
)

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
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "woadult_scv.pl.velocity_embedding_stream_scpred_prediction.svg",
    bbox_inches="tight",
    backend = "cairo",
    transparent=True,
)

################################################################################################


adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG_ldata_dynamics.h5ad"
)
adata.obs["Days"] = adata.obs["Days"].astype(float)
adata.obs["Days"] = adata.obs["Days"].replace(-1, np.nan)

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)
adata.obs["scpred_prediction"] = pd.Categorical(
    list(adata.obs["scpred_prediction"]),
    categories=[
        "RPC",
        "RGC",
        "Cone",
        "HC",
        "AC",
        "Rod",
        "BC",
        "MG",
    ],
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    legend_fontsize=20,
    title="",
    color="Days",
    size=10,
    colorbar=True,
    dpi=600,
    figsize=(10, 10),
    linewidth=5,
    alpha=1,
    legend_loc="right margin",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "wadult_scv.pl.velocity_embedding_stream_Days.svg",
    bbox_inches="tight",
    backend = "cairo",
    transparent=True,
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    legend_fontsize=20,
    title="",
    color="scpred_prediction",
    size=10,
    colorbar=True,
    dpi=600,
    figsize=(10, 10),
    linewidth=5,
    alpha=1,
    legend_loc="right margin",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "wadult_scv.pl.velocity_embedding_stream_scpred_prediction.svg",
    bbox_inches="tight",
    backend = "cairo",
    transparent=True,
)
