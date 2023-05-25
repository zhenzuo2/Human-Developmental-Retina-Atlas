import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/HC.h5ad"
)
meta_file = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/HC_subtype.csv"
)
meta_file.index = meta_file["Unnamed: 0"]
common_cells = [x for x in meta_file.index if x in adata_result.obs.index]
meta_file = meta_file.loc[common_cells, :]
adata_result.obs["subclass"] = adata_result.obs["majorclass"].astype(str)
adata_result.obs["majorclass"] = adata_result.obs["majorclass"].astype(str)

for sublcass in set(meta_file.subclass):
    adata_result.obs.loc[
        meta_file.loc[meta_file.subclass == sublcass, :].index, "subclass"
    ] = sublcass
for majorclass in set(meta_file.majorclass):
    adata_result.obs.loc[
        meta_file.loc[meta_file.majorclass == majorclass, :].index, "majorclass"
    ] = majorclass

adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
adata_result = adata_result[adata_result.obs.majorclass!="HC"]

mv.velocity_embedding_stream(adata_result, basis="umap")

mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="majorclass",
    alpha=1,
    size=50,
    legend_fontsize=20,
    legend_loc="on data",
    frameon=False,
    legend_fontoutline=True,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/HC_mv.velocity_embedding_stream.majorclass.svg",
    dpi=600,
)

mv.velocity_embedding_stream(
    adata_result, basis="umap", color="Days", size=50, alpha=1, frameon=False,legend_fontoutline=True,legend_fontsize=20,cmap = "inferno"
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/HC_mv.velocity_embedding_stream.Days.svg",
    dpi=600,
)

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)
mv.velocity_embedding_stream(
    adata_result, basis="umap", color="latent_time", size=50, alpha=1, frameon=False,legend_fontoutline=True,legend_fontsize=20,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/HC_mv.velocity_embedding_stream.latent_time.svg",
    dpi=600,
)



