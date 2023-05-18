# Import packages
import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

scv.set_figure_params(dpi=600, dpi_save=600)

gs = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/gene_score.h5ad"
)
adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad"
)
scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)

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
}

df = pd.DataFrame(gs.obs)
df.loc[:, "cell_id"] = list(df.index.values)
grouped_data = df.groupby("majorclass")
# Define the number of rows to downsample
num_rows = 10000
# Sample the same number of rows from each group
downsampled_data = grouped_data.apply(
    lambda x: x.sample(n=num_rows, random_state=42, replace=True)
)
# Reset the index
downsampled_data = downsampled_data.reset_index(drop=True)

sc.set_figure_params(scanpy=True, dpi_save=600, fontsize=25)
gs_subset = gs[downsampled_data.cell_id]
gs_subset.obs["majorclass"] = pd.Categorical(
    list(gs_subset.obs["majorclass"]),
    categories=["BC","Cone","Rod","MG","RPC","RGC","AC","HC"],
)
sc.pl.heatmap(
    gs_subset,
    Markers,
    groupby="majorclass",
    cmap="viridis",
    dendrogram=False,
    show_gene_labels=True,
    vmax=np.quantile(gs_subset.X, 0.95),)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/gene_score_heatmap.svg",
    bbox_inches="tight",
    transparent=True,
)
##########################################################################################################################
adata_subset = adata[downsampled_data.cell_id]
ax = sc.pl.heatmap(
    adata_subset,
    Markers,
    groupby="majorclass",
    dendrogram=True,
    show_gene_labels=True,
    vmax=np.quantile(adata_subset.X.toarray(), 0.99),
    cmap = "plasma"
)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/gene_expression_heatmap.svg",
    bbox_inches="tight",
    transparent=True,
)
