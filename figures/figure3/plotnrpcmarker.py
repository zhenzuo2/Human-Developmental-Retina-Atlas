# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns

sc.set_figure_params(transparent=True, fontsize=80)
NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC = NRPC[NRPC.obs.leiden!="12"]
NRPC.obs["subclass"] = np.nan

clusters = ["2", "13"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["0", "1", "18"]
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
    NRPC,
    size=40,
    color="subclass",
    title="",
    frameon=False,
    legend_loc="right margin",
    palette={
        "AC": "#1f77b4",
        "BC": "#ff7f0e",
        "Cone": "#2ca02c",
        "HC": "#d62728",
        "RGC": "#7f7f7f",
        "Rod": "#bcbd22",
    },
),
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass.svg",
    dpi=600,
    bbox_inches="tight",
    backend = "cairo"
)
plt.clf()

sc.pp.normalize_total(NRPC)
sc.pp.log1p(NRPC)

NRPC.obs["Days"] = NRPC.obs["Days"].astype(float)
NRPC.obs["PCW"] = NRPC.obs["Days"] / 7

sc.pl.umap(NRPC, size=40, color="Days", title="", frameon=False, cmap="cividis"),
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_Days.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

NRPC.obs["Weeks"] = NRPC.obs.Days.map(
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
        palette={
            "AC": "#1f77b4",
            "BC": "#ff7f0e",
            "Cone": "#2ca02c",
            "HC": "#d62728",
            "RGC": "#7f7f7f",
            "Rod": "#bcbd22",
        },
    ),
    fig = plt.gcf()
    fig.set_size_inches(10, 10)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_" + sp + ".svg",
        dpi=600,
        bbox_inches="tight",
    )
    plt.clf()

cells = []
for sp in [
    "AC",
    "BC",
    "Cone",
    "HC",
    "RGC",
    "Rod",
]:
    cells = cells + [NRPC[NRPC.obs.subclass == sp].obs.index[0]]


for region in set(NRPC.obs.Region):
    for weeks in set(NRPC.obs.Weeks):
        NRPC.obs["temp"] = np.nan
        subset = (NRPC.obs.Region == region) & (NRPC.obs.Weeks == weeks)
        NRPC.obs.loc[cells, "temp"] = NRPC.obs.loc[cells, "subclass"]
        NRPC.obs.loc[subset, "temp"] = NRPC.obs.loc[subset, "subclass"]
        NRPC.obs.loc[subset, "temp"] = NRPC.obs.loc[subset, "temp"].astype(str)
        print(NRPC.obs.loc[subset, "temp"])
        sc.pl.umap(
            NRPC,
            color="temp",
            size=80,
            title="",
            frameon=False,
            legend_loc="None",
            na_color="white",
            palette={
                "AC": "#1f77b4",
                "BC": "#ff7f0e",
                "Cone": "#2ca02c",
                "HC": "#d62728",
                "RGC": "#7f7f7f",
                "Rod": "#bcbd22",
                "nan": "lightgray",
            },
        )
        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        plt.savefig(
            "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/overall_umap_by_cell_type_"
            + region
            + "_"
            + weeks
            + ".svg",
            dpi=600,
            bbox_inches="tight",
            transparent=True,
        )


sc.pp.normalize_total(NRPC)
sc.pp.log1p(NRPC)


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
plt.rcParams.update(
    {
        "font.size": 40,
    }
)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_DEGs_dot_plot.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# AC
sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["NEUROD4", "PAX6", "PTF1A", "PRDM13"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_AC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# RGC
sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["POU4F2", "ISL1", "ATOH7", "POU4F1"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_RGC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# HC
sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["ONECUT1", "ONECUT2", "ONECUT3", "PROX1"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_HC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# BC
sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["OTX2", "VSX2", "PRDM1", "VSX1"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_BC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# Rod
sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["NRL", "CRX", "NR2E3", "OTX2"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_Rod.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

# Cone
sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["THRB", "CRX", "PRDM1", "OTX2"],
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(40, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_Cone.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["AC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/AC_fate_p.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["BC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/BC_fate_p.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["Rod"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/Rod_fate_p.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["Cone"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/Cone_fate_p.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["RGC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/RGC_fate_p.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["HC"],
    cmap="plasma",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/HC_fate_p.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()
