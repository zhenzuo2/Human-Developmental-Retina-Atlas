# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import pandas as pd

sc.set_figure_params(transparent=True, fontsize=80)
NRPC = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.h5ad"
)
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/to_terminal_states.csv"
)
df.index = df["Unnamed: 0"].values

for x in ["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]:
    NRPC.obs[x] = df.loc[NRPC.obs.index, x]

# AC
sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["NEUROD4", "PAX6", "PTF1A", "PRDM13"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(45, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_subclass_AC.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

# RGC
sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["POU4F2", "ISL1", "ATOH7", "POU4F1"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(45, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_subclass_RGC.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

# HC
sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["ONECUT1", "ONECUT2", "ONECUT3", "PROX1"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(45, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_subclass_HC.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

# BC
sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["OTX2", "VSX2", "PRDM1", "VSX1"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(45, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_subclass_BC.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

# Rod
sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["NRL", "CRX", "NR2E3", "OTX2"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(45, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_subclass_Rod.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

# Cone
sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["THRB", "CRX", "PRDM1", "OTX2"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(45, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_subclass_Cone.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["AC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/AC_fate_p.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["BC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/BC_fate_p.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["Rod"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/Rod_fate_p.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["Cone"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/Cone_fate_p.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["RGC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/RGC_fate_p.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    #legend_loc=None,
    #colorbar_loc=None,
    color=["HC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/HC_fate_p.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()


for sp in [
    "AC",
    "BC",
    "Cone",
    "HC",
    "RGC",
    "Rod",
]:
    NRPC.obs["temp"] = np.nan
    subset = NRPC.obs.index[NRPC.obs.subclass == sp]
    NRPC.obs.loc[subset, "temp"] = sp
    sc.pl.umap(
        NRPC,
        size=40,
        color="temp",
        title="",
        frameon=False,
        legend_loc=None,
        palette={
            "MG": "#9467bd",
            "Rod": "#17becf",
            "BC": "#bcbd22",
            "RGC": "#d62728",
            "NRPC": "#d3d3d3",
            "Cone": "#e377c2",
            "HC": "#2ca02c",
            "PRPC": "#1f77b4",
            "AC": "#8c564b",
        },
    ),
    fig = plt.gcf()
    fig.set_size_inches(10, 10)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_" + sp + ".tiff",
        dpi=300,
        bbox_inches="tight",
    )
    plt.clf()
