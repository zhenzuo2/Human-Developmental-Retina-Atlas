# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

plt.rcParams.update({"font.size": 20})
plt.rcParams["font.family"] = "Arial"
sc.pl.violin(
    adata,
    groupby="majorclass",
    keys="OTX2",
    rotation=90,
    order=["Rod", "Cone", "BC", "NRPC", "PRPC", "AC", "MG", "HC", "RGC"],
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
    },
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.xlabel("")
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/sc.pl.violin_OTX2.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
