# Import packages
import scanpy as sc
import scvi
import os
import scvelo as scv
import anndata
import sys
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
import numpy as np
import tempfile

scv.settings.verbosity = 3

scv.set_figure_params(dpi=200, dpi_save=600)

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad"
)

Markers = {
    "BC": [
        "VSX1",
        "VSX2",
        "OTX2",
        "GRM6",
        "PRKCA",
        "LHX4",
        "PROX1",
        "PCP4",
        "PCP2",
        "TRPM1",
        "PRDM8",
    ],
    "Cone": [
        "PROM1",
        "CRX",
        "ARR3",
        "GNAT2",
        "THRB",
        "OPN1SW",
        "PDE6H",
        "GADD45G",
        "NEUROD1",
        "RXRG",
        "DCT",
        "PRDM1",
        "RAX2",
        "CRABP2",
        "HOTAIRM1",
    ],
    "Rod": [
        "PROM1",
        "CRX",
        "RCVRN",
        "OTX2",
        "RHO",
        "NR2E3",
        "GNAT1",
        "NRL",
        "GADD45G",
        "NEUROD1",
        "RXRG",
        "DCT",
        "PRDM1",
        "RAX2",
        "CRABP2",
        "PRPH2",
        "SAG",
    ],
    "RGC": [
        "POU4F2",
        "RBPMS",
        "NEFM",
        "GAP43",
        "POU4F1",
        "ELAVL4",
        "POU6F2",
        "ISL1",
        "NHLH2",
        "RXRG",
        "EBF1",
        "EBF3",
        "MYC",
    ],
    "AC": [
        "SLC6A9",
        "GAD1",
        "SLC32A1",
        "TFAP2B",
        "GAD2",
        "SLC18A3",
        "LHX9",
        "MEIS2",
        "TFAP2C",
        "TFAP2A",
    ],
    "HC": [
        "ONECUT1",
        "ONECUT2",
        "ONECUT3",
        "TFAP2B",
        "LHX1",
        "TFAP2A",
        "ESRRB",
    ],
    "MG": [
        "SLC1A3",
        "SLN",
        "RLBP1",
        "SOX2",
        "NFIA",
        "CRYM",
        "CLU",
        "LINC00461",
    ],
    "RPC": [
        "VIM",
        "SOX2",
        "SFRP2",
        "MKI67",
        "UBE2C",
        "FGF19",
        "CCND1",
        "ID3",
    ],
}

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
scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)
df = adata.obs
df["cell_id"] = df.index.values
grouped_data = df.groupby("scpred_prediction")
# Define the number of rows to downsample
num_rows = 10000
# Sample the same number of rows from each group
downsampled_data = grouped_data.apply(
    lambda x: x.sample(n=num_rows, random_state=42, replace=True)
)
# Reset the index
downsampled_data = downsampled_data.reset_index(drop=True)

sc.set_figure_params(scanpy=True, dpi_save=600, fontsize=25)
ax = sc.pl.heatmap(
    adata[downsampled_data.cell_id],
    Markers,
    groupby="scpred_prediction",
    cmap="viridis",
    dendrogram=True,
    show_gene_labels=True,
    save="heatmap.svg",
)
sc.pl.umap(
    adata,
    color="CCND1",
    size=1,
    legend_loc="right margin",
    title="CCND1",
    vmax=1,
    save="CCND1.svg",
    frameon=False,
)
sc.pl.umap(
    adata,
    color="ID3",
    size=1,
    legend_loc="right margin",
    title="ID3",
    vmax=1,
    save="ID3.svg",
    frameon=False,
)
sc.pl.umap(
    adata,
    color="NFIA",
    size=1,
    legend_loc="right margin",
    title="NFIA",
    vmax=3,
    save="NFIA.svg",
    frameon=False,
)
sc.pl.umap(
    adata,
    color="ASCL1",
    size=1,
    legend_loc="right margin",
    title="ASCL1",
    vmax=1,
    save="ASCL1.svg",
    frameon=False,
)