import os
import numpy as np
import pandas as pd
import scanpy as sc
import multivelo as mv
import sys
import matplotlib.pyplot as plt

adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC_MG_NRPC/adata_umap.h5ad"
)
sc.pl.umap(
adata,
color="majorclass",
size=80,
title="",
frameon=False,
legend_loc="None",
na_color="lightgray",
palette={
    "PRPC": "#1f77b4",
    "MG": "#ff7f0e",
    "NRPC":"#2ca02c",
},
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
"/storage/singlecell/zz4/fetal_snakemake/figures/figure6/overall_umap_by_majorclass.png",
dpi=600,
bbox_inches="tight",
transparent=True,
)

adata.obs["Weeks"] = adata.obs.Days.map(
    {
        70: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW16",
        103: "PCW16",
        116: "PCW16",
        120: "PCW16",
        136: "PCW20",
        137: "PCW20",
        141: "PCW20",
        142: "PCW20",
        162: "PCW23",
        165: "PCW23",
    }
)

for region in set(adata.obs.Region):
    for weeks in set(adata.obs.Weeks):
        adata.obs["temp"] = np.nan
        subset = (adata.obs.Region == region) & (adata.obs.Weeks == weeks)
        adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "majorclass"]
        adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "temp"].astype(str)
        print(adata.obs.loc[subset, "temp"])
        sc.pl.umap(
            adata,
            color="temp",
            size=80,
            title="",
            frameon=False,
            legend_loc="None",
            na_color="lightgray",
            palette={
                "PRPC": "#1f77b4",
                "MG": "#ff7f0e",
                "NRPC":"#2ca02c",
            },
        )
        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        plt.savefig(
            "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/overall_umap_by_majorclass_"
            + region
            + "_"
            + weeks
            + ".png",
            dpi=600,
            bbox_inches="tight",
            transparent=True,
        )
