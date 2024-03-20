# Import packages
import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import anndata
from scipy.stats.stats import pearsonr
from csaps import csaps
import matplotlib as mpl

mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["scatter.edgecolors"] = "none"
font = {"family": "Arial", "size": 25}
plt.rc("font", **font)

dev = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
dev.obs["Region"] = dev.obs["Region"].replace({"Peripheral":"Periphery"})
adult = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/all/major_clean_major_scvi_Cluster_clean.h5ad"
)
adult.obs["sampleid"] = "Adult"
adult.obs["Region"] = "Adult"
adult.obs["Days"] = "Adult"

adata = anndata.concat([adult, dev])

adata.obs["Weeks"] = adata.obs.Days.map(
    {
        59: "PCW8",
        70: "PCW10",
        76: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW15",
        103: "PCW15",
        116: "PCW15",
        137: "PCW19",
        141: "PCW19",
        142: "PCW19",
        162: "PCW23",
        165: "PCW23",
        "Adult": "Adult",
    }
)

adata = adata[
    ~(
        (adata.obs.Days == 142)
        & (adata.obs.Region == "Macula")
        & (adata.obs.majorclass.isin(["RGC", "HC"]))
    )
]

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=True)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.tl.pca(adata)


def get_pca(adata):
    return np.mean(adata.obsm["X_pca"], axis=0)


def cal_cor(adata, sam1, sam2, region):
    if (adata[(adata.obs.Days == sam1) & (adata.obs.Region == region)].shape[0]) > 0:
        pc1 = get_pca(adata[(adata.obs.Days == sam1) & (adata.obs.Region == region)])
        pc2 = get_pca(adata[(adata.obs.Days == sam2)])
        return pearsonr(pc1, pc2)[0]


def smooth(df):
    temp = df[df.Region == "Periphery"]
    x = np.array(temp.Days)
    y = np.array(temp.Cor)
    xs = np.linspace(min(x), max(x), len(x))
    ys = csaps(x, y, xs)
    df.loc[df.Region == "Periphery", "ys"] = ys[0]

    temp = df[df.Region == "Macula"]
    x = np.array(temp.Days)
    y = np.array(temp.Cor)
    xs = np.linspace(min(x), max(x), len(x))
    ys = csaps(x, y, xs)
    df.loc[df.Region == "Macula", "ys"] = ys[0]
    return df


def plot_curve(adata, cell_type):
    temp = adata[adata.obs.majorclass == cell_type]
    df = pd.DataFrame(columns=["Days", "Region", "Cor"])
    for x in [
        59,
        70,
        76,
        79,
        87,
        91,
        100,
        103,
        116,
        137,
        141,
        142,
        162,
        165,
    ]:
        df = df.append(
            {"Days": x, "Region": "Macula", "Cor": cal_cor(temp, x, "Adult", "Macula")},
            ignore_index=True,
        )
        df = df.append(
            {
                "Days": x,
                "Region": "Periphery",
                "Cor": cal_cor(temp, x, "Adult", "Periphery"),
            },
            ignore_index=True,
        )
    df = df.dropna()
    df = smooth(df)
    df.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/Cell Maturation.csv")
    if cell_type == "RGC":
        condition1 = (df["Region"] == "Macula") & (df["Days"] >= 70)
        condition2 = (df["Region"] == "Periphery") & (df["Days"] >= 70)
        df = df[condition1 | condition2]
    if cell_type == "HC":
        condition1 = (df["Region"] == "Macula") & (df["Days"] >= 76)
        condition2 = (df["Region"] == "Periphery") & (df["Days"] >= 79)
        df = df[condition1 | condition2]
    if cell_type == "Cone":
        condition1 = (df["Region"] == "Macula") & (df["Days"] >= 76)
        condition2 = (df["Region"] == "Periphery") & (df["Days"] >= 91)
        df = df[condition1 | condition2]
    if cell_type == "AC":
        condition1 = (df["Region"] == "Macula") & (df["Days"] >= 100)
        condition2 = (df["Region"] == "Periphery") & (df["Days"] >= 103)
        df = df[condition1 | condition2]
    if cell_type == "Rod":
        condition1 = (df["Region"] == "Macula") & (df["Days"] >= 103)
        condition2 = (df["Region"] == "Periphery") & (df["Days"] >= 137)
        df = df[condition1 | condition2]
    if cell_type == "BC":
        condition1 = (df["Region"] == "Macula") & (df["Days"] >= 103)
        condition2 = (df["Region"] == "Periphery") & (df["Days"] >= 142)
        df = df[condition1 | condition2]
    if cell_type == "MG":
        condition1 = (df["Region"] == "Macula") & (df["Days"] >= 76)
        condition2 = (df["Region"] == "Periphery") & (df["Days"] >= 142)
        df = df[condition1 | condition2]

    fig, ax = plt.subplots()
    for color, group in df.groupby("Region"):
        if color == "Macula":
            ax.plot(
                group["Days"],
                group["ys"],
                marker="o",
                linestyle="-",
                label=color,
                color="#F8766D",
            )
        if color == "Periphery":
            ax.plot(
                group["Days"],
                group["ys"],
                marker="o",
                linestyle="-",
                label=color,
                color="#00BFC4",
            )
    # Adding labels and title
    ax.set_xlabel("Postconceptional age (Days)")
    ax.set_ylabel("Cell Maturation Index")
    ax.set_ylim([0, 1])
    ax.set_xlim([59, 175])
    ax.set_title(cell_type, fontsize=25)
    # Adding legend
    ax.legend()
    # Display the plot
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/Cell Maturation Index_"
        + cell_type
        + ".tiff",
        dpi=300,
        bbox_inches="tight",
        transparent=True,
    )

for cell_type in ["AC", "BC", "Cone", "HC", "RGC", "Rod","MG"]:
    plot_curve(adata, cell_type)
