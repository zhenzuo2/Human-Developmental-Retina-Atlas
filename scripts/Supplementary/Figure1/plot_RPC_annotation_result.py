import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


adata_result = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

matplotlib.rcParams.update({"font.size": 12})
sc.pl.umap(adata_result, color="majorclass", size=10, title="")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/majorclass_annotation.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)

matplotlib.rcParams.update({"font.size": 6})
adata_result.obs["majorclass"] = adata_result.obs.majorclass.astype(str)
adata_result.obs["leiden"] = adata_result.obs.leiden.astype(str)
adata_result.obs.loc[
    ~adata_result.obs.majorclass.isin(["PRPC", "NRPC", "MG"]), "majorclass"
] = "Non-Progenitor Cells"
adata_result.obs.loc[
    ~adata_result.obs.majorclass.isin(["PRPC", "NRPC", "MG"]), "leiden"
] = "Non-Progenitor Cells"
sc.pl.umap(adata_result, color="leiden", size=10, title="", legend_loc="on data")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/RPC_annotation_leiden.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)

matplotlib.rcParams.update({"font.size": 12})
sc.pl.umap(adata_result, color="majorclass", size=10, title="",palette={
                "MG": "#9467bd",
                "NRPC": "#ff7f0e",
                "PRPC": "#1f77b4",
                "Non-Progenitor Cells": "lightgray"
            },frameon=False,)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/RPC_annotation.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)
