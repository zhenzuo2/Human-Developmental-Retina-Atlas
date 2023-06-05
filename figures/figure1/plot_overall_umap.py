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
    backend = "cairo"
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
    backend = "cairo"
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
    backend = "cairo"
)
#############################################################################################################################
## Plot subtype
adata.obs["subclass"] = ""
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
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_subtype.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend = "cairo"
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
            + ".svg",
            dpi=600,
            bbox_inches="tight",
            transparent=True,
            backend = "cairo"
        )
