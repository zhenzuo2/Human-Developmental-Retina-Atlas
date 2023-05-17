# Import packages
# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns

NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC.obs["subclass"] = np.nan

clusters = ["2", "13"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["1"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["9", "10"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["4"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["17"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["7", "14"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"

sc.pl.umap(
    NRPC, size=40, color="subclass", title="", frameon=False, legend_loc="right margin"
),
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pp.normalize_total(NRPC)
sc.pp.log1p(NRPC)

NRPC.obs["Days"] = NRPC.obs["Days"].astype(float)
NRPC.obs["PCW"] = NRPC.obs["Days"] / 7
# Create some sample data

# Create a boxplot using Seaborn
sns.boxplot(
    data=NRPC.obs,
    x="subclass",
    y="PCW",
    hue="Region",
    order=["RGC", "Cone", "HC", "AC", "Rod", "BC"],
    showfliers=False,
)
plt.title("")
plt.xlabel("Precursor Group")
plt.ylabel("PCW")

# Move the legend outside the figure
plt.legend(title="Region", loc="best", bbox_to_anchor=(1, 0.5))
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_Days.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()


sc.tl.rank_genes_groups(NRPC, "subclass")
sc.pl.rank_genes_groups_dotplot(NRPC, n_genes=10)
fig = plt.gcf()
fig.set_size_inches(50, 10)
plt.rcParams.update({'font.size': 30,})
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_DEGs_dot_plot.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# AC
scv.pl.umap(NRPC, size=40, color=["NEUROD4", "PAX6", "PTF1A", "PRDM13"])
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_AC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# RGC
scv.pl.umap(NRPC, size=40, color=["POU4F2", "ISL1", "ATOH7", "POU4F1"])
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_RGC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# HC
scv.pl.umap(NRPC, size=40, color=["ONECUT1", "ONECUT2", "ONECUT3", "PROX1"])
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_HC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# BC
scv.pl.umap(NRPC, size=40, color=["OTX2", "VSX2", "PRDM1", "VSX1"])
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_BC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# Rod
scv.pl.umap(NRPC, size=40, color=["NRL", "CRX", "NR2E3", "OTX2"])
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_Rod.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# Cone
scv.pl.umap(NRPC, size=40, color=["THRB", "CRX", "PRDM1", "OTX2"])
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_Cone.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()
