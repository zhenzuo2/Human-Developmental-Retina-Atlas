# Import packages
import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
scv.set_figure_params(dpi=600, dpi_save=600)

gs = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/gene_score.h5ad"
)
adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)

gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]

cells = [x for x in adata.obs.index if x in gs.obs.index]
adata = adata[cells]
gs = gs[cells]
adata.obs = pd.DataFrame(adata.obs)
adata.obs["majorclass"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)
gs.obs = adata.obs
gs.obsm['X_umap'] = adata.obsm['X_umap']

sc.pp.normalize_per_cell(gs)
sc.pp.log1p(gs)

sc.pl.umap(
    gs[gs.obs.Region =="Macula"],
    vmax = 1.5,
    color="CNGB3",
    size=10,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/CNGB3_Macula_gs.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)

sc.pl.umap(
    gs[gs.obs.Region =="Peripheral"],
    vmax = 1.5,
    color="CNGB3",
    size=10,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/CNGB3_Peripheral_gs.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)