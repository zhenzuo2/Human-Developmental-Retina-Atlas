# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
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
#############################################################################################################################
## Plot subtype
adata.obs["subclass"] = ""
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs.loc[adata.obs.majorclass == "PRPC", "subclass"] = "PRPC"
adata.obs.loc[adata.obs.majorclass == "NRPC", "subclass"] = "NRPC"
adata.obs.loc[adata.obs.majorclass == "MG", "subclass"] = "MG"
for file in [
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/AC_subtype.csv",
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/BC_subtype.csv",
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/Cone_subtype.csv",
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/HC_subtype.csv",
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/RGC_subtype.csv",
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/Rod_subtype.csv",
]:
    df = pd.read_csv(file)
    df.index = df["Unnamed: 0"].values
    adata.obs.loc[df.index, "subclass"] = df.subclass
    adata.obs.loc[df.index, "majorclass"] = df.majorclass

adata = adata[adata.obs.subclass != "", :]

scv.pl.umap(
    adata,
    color="scpred_prediction",
    size=5,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_cell_type.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)

scv.pl.umap(
    adata,
    color="Region",
    size=5,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_region.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)

adata.obs["Days"] = adata.obs["Days"].astype(float)
scv.pl.umap(
    adata,
    color="Days",
    size=5,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_days.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)

sc.pl.umap(
    adata,
    color="subclass",
    size=5,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
    palette = ['#88ccee', '#cc6677', '#ddcc77', '#117733', '#332288', '#aa4499',
             '#44aa99', '#999933', '#882255', '#661100', '#888888']
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_subtype.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)

plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "Arial"
sc.pl.umap(
    adata,
    color="majorclass",
    size=5,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
    palette = 'tab20'
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_majortype.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
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
cells = []
for sp in set(adata.obs.scpred_prediction):
    cells = cells + [adata[adata.obs.scpred_prediction == sp].obs.index[0]]

for region in set(adata.obs.Region):
    for weeks in set(adata.obs.Weeks):
        adata.obs["temp"] = np.nan
        subset = (adata.obs.Region == region) & (adata.obs.Weeks == weeks)
        adata.obs.loc[cells, "temp"] = adata.obs.loc[cells, "scpred_prediction"]
        adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "scpred_prediction"]
        adata.obs["temp"] = pd.Categorical(
            list(adata.obs["temp"]),
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
        sc.pl.umap(
            adata,
            color="temp",
            size=10,
            title="",
            frameon=False,
            legend_loc="None",
            palette={
                "RPC": "#1f77b4",
                "RGC": "#ff7f0e",
                "Cone": "#2ca02c",
                "HC": "#d62728",
                "AC": "#9467bd",
                "Rod": "#8c564b",
                "BC": "#e377c2",
                "MG": "#7f7f7f",
            },
        )
        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        plt.savefig(
            "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_cell_type_"
            + region
            + "_"
            + weeks
            + ".png",
            dpi=600,
            bbox_inches="tight",
            transparent=True,
            #backend="cairo",
        )


################################################################################################################################################################
################################################################################################################################################################
week_set = []
for weeks in ['PCW10', 'PCW13','PCW16', 'PCW20', 'PCW23']:
    week_set =  [weeks]
    adata.obs["temp"] = ""
    subset = adata.obs.Weeks.isin(week_set)
    adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "scpred_prediction"]
    adata.obs["temp"] = pd.Categorical(
        list(adata.obs["temp"]),
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
    sc.pl.umap(
        adata[(adata.obs.Weeks.isin(week_set))&(adata.obs.Region =="Macula"),],
        color="temp",
        size=20,
        title="",
        frameon=True,
        legend_loc="None",
        palette={
            "RPC": "#1f77b4",
            "RGC": "#ff7f0e",
            "Cone": "#2ca02c",
            "HC": "#d62728",
            "AC": "#9467bd",
            "Rod": "#8c564b",
            "BC": "#e377c2",
            "MG": "#7f7f7f",
        },
        outline_color = "white",
    )
    fig = plt.gcf()
    fig.set_size_inches(10, 10)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_cell_type_"
        + "_"
        + weeks
        + ".png",
        dpi=600,
        bbox_inches="tight",
        transparent=True,
        #backend="cairo",
    )

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "Arial"
sc.pl.violin(
    adata,
    groupby='majorclass',
    keys="OTX2",
    title="",
    rotation=90,
    order = ["Rod", "Cone","BC", "NRPC", "PRPC", "AC","MG", "HC", "RGC"]
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.xlabel('')
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_OTX2.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="GRB10",
    size=5,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/GRB10.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)