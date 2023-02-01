# Import packages
import scanpy as sc
import scvi
import os
import scvelo as scv
import anndata
import sys
import pandas as pd
import matplotlib.pyplot as plt
import cellrank as cr
import seaborn as sns
import numpy as np
import tempfile

# import CellRank kernels and estimators
from cellrank.external.kernels import WOTKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA

# set verbosity levels
cr.settings.verbosity = 2
scv.settings.verbosity = 3

scv.set_figure_params(dpi=200, dpi_save=600)

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000.h5ad"
)
adata.obsm["X_umap"]
umap = pd.DataFrame(adata.obsm["X_umap"], columns=["x", "y"])
umap.index = adata.obs.index
umap

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad"
)

adata.obsm["X_umap"] = umap.loc[
    adata.obs.index,
].values

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "AC Precursor": "AC",
        "BC Precursor": "BC",
        "Cone Precursor": "Cone",
        "GABAergic": "AC",
        "Glycinergic": "AC",
        "HC0": "HC",
        "HC1": "HC",
        "MG": "MG",
        "ML_Cone": "Cone",
        "NRPC": "RPC",
        "OFF-BC": "BC",
        "OFF_MGC": "RGC",
        "ON-BC": "BC",
        "ON_MGC": "RGC",
        "PRPC": "RPC",
        "RBC": "BC",
        "RGC Precursor": "RGC",
        "Rod": "Rod",
        "Rod Precursor": "Rod",
        "SACs": "AC",
        "S_Cone": "Cone",
        "dual ACs": "AC",
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
scv.pl.umap(
    adata,
    color="scpred_prediction",
    size=1,
    legend_loc="right margin",
    title="",
    save="figure1b.svg",
)

scv.set_figure_params(dpi=200, dpi_save=600, fontsize=5)
scv.pl.umap(
    adata,
    color="subclass",
    size=1,
    legend_loc="on data",
    title="",
    save="figure1d.svg",
    legend_fontweight="normal",
)
scv.set_figure_params(fontsize=14)

scv.pl.umap(
    adata[adata.obs.Region.isin(["Macula", "Peripheral"])],
    color="Region",
    size=1,
    legend_loc="right margin",
    title="",
    save="figure1e.svg",
)

adata.obs["Days"] = adata.obs["Days"].astype(float)
scv.pl.umap(
    adata[adata.obs.Region.isin(["Macula"])],
    color="Days",
    size=1,
    legend_loc="right margin",
    title="",
    save="figure1f.svg",
    color_map="plasma",
)


adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_wadult_umap_10000.h5ad"
)
adata = adata[~adata.obs.sampleid.isin(["17w1d_I_Ret"])]
adata.obs["batch"] = adata.obs["batch"].replace(
    {"Query": "Fetal", "Reference": "Adult"}
)

scv.pl.umap(
    adata,
    color="batch",
    size=1,
    legend_loc="right margin",
    title="",
    save="figure1g.svg",
)
adata.obs["majorclass"] = adata.obs["majorclass"].replace({"MG": "RPC/MG"})
scv.pl.umap(
    adata,
    color="majorclass",
    size=1,
    # legend_loc="right margin",
    title="",
    save="figure1f.svg",
)
scv.pl.umap(
    adata[
        adata.obs.batch == "Fetal",
    ],
    color="Data Type",
    size=0.5,
    legend_loc="right margin",
    title="",
    sort_order=False,
    save="figure1h.svg",
)