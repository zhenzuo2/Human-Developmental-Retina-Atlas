# Import packages
import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import magic


def smooth(adata):
    magic_op = magic.MAGIC()
    magic_op.set_params(n_jobs=10)
    emt_magic = magic_op.fit_transform(adata.X, genes="all_genes")
    emt_magic = magic_op.transform(genes="all_genes")
    adata.X = np.array(emt_magic)
    return adata

scv.set_figure_params(dpi=300, dpi_save=300)

gs = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score.h5ad"
)
adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

sc.pp.normalize_total(gs, target_sum=1e6)
sc.pp.log1p(gs)

gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]

cells = [x for x in adata.obs.index if x in gs.obs.index]
adata = adata[cells]
gs = gs[cells]
adata.obs = pd.DataFrame(adata.obs)

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
    "PRPC": [
        "VIM",
        "SOX2",
        "SFRP2",
        "MKI67",
        "UBE2C",
        "FGF19",
        "CCND1",
        "ID3",
        "SFRP2",
        "NFIA",
        "DLX1",
        "ASCL1",
        "DLX2",
        "OTX2",
        "ONECUT2",
        "SOX4",
        "ATOH7",
    ],
    "NRPC": [
        "DLX1",
        "ASCL1",
        "DLX2",
        "OTX2",
        "ONECUT2",
        "SOX4",
        "ATOH7",
        "OLIG2",
        "NEUROG2",
        "BTG2",
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
num_rows = 2000
# Sample the same number of rows from each group
downsampled_data = grouped_data.apply(
    lambda x: x.sample(n=num_rows, random_state=42, replace=True)
)
# Reset the index
downsampled_data = downsampled_data.reset_index(drop=True)

sc.set_figure_params(scanpy=True, dpi_save=600, fontsize=25)
gs_subset = gs[downsampled_data.cell_id]
all_values = [value for values in Markers.values() for value in values]
all_values = list(set(all_values))
gs_subset = gs_subset[:, all_values]
sc.pp.scale(gs_subset)
gs_subset = smooth(gs_subset)
gs_subset.X = np.array(gs_subset.X)
gs_subset.obs["majorclass"] = pd.Categorical(
    list(gs_subset.obs["majorclass"]),
    categories=["BC", "Cone", "Rod", "MG", "PRPC", "NRPC", "RGC", "AC", "HC"],
)
gs_subset.write(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/gene_score_heatmap_markers.h5ad"
)
sc.pl.heatmap(
    gs_subset,
    Markers,
    groupby="majorclass",
    cmap="plasma",
    dendrogram=False,
    show_gene_labels=True,
    standard_scale="obs",
    vmax = 0.7,
    vmin = 0.4
)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/gene_score_heatmap_markers.tiff",
    bbox_inches="tight",
    transparent=True,
)
##########################################################################################################################
adata_subset = adata[downsampled_data.cell_id, all_values]
sc.pp.scale(adata_subset)
adata_subset = smooth(adata_subset)
adata_subset.X = np.array(adata_subset.X)
adata_subset.obs["majorclass"] = pd.Categorical(
    list(adata_subset.obs["majorclass"]),
    categories=["BC", "Cone", "Rod", "MG", "PRPC", "NRPC", "RGC", "AC", "HC"],
)
adata_subset.write(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/gene_expression_heatmap_markers.h5ad"
)
sc.pl.heatmap(
    adata_subset,
    Markers,
    groupby="majorclass",
    dendrogram=False,
    show_gene_labels=True,
    cmap="viridis",
    standard_scale="obs",
    vmax = 0.5,
    vmin = 0.1
)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/gene_expression_heatmap_markers.tiff",
    bbox_inches="tight",
    transparent=True,
)
