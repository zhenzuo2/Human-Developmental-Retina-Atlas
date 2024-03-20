import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["scatter.edgecolors"] = "none"
font = {"family": "Arial", "size": 12}
plt.rc("font", **font)
plt.figure(figsize=(5, 5))


def plot_umap(x, y, c, size, save):
    plt.scatter(list(x), list(y), c=c, s=size, alpha=1)

    plt.title("")
    plt.xlabel("")
    plt.ylabel("")
    plt.box(False)
    plt.tick_params(
        top="off",
        bottom="off",
        left="off",
        right="off",
        labelleft="off",
        labelbottom=",",
    )
    plt.xticks([])
    plt.yticks([])

    plt.savefig(
        save,
        dpi=300,
        transparent=True,
    )


df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/ATAC_umap/ATAC_UMAPHarmony_cor.csv"
)
pairwise_dict = {
    "MG": "#9467bd",
    "Rod": "#17becf",
    "BC": "#bcbd22",
    "RGC": "#d62728",
    "NRPC": "#ff7f0e",
    "Cone": "#e377c2",
    "HC": "#2ca02c",
    "PRPC": "#1f77b4",
    "AC": "#8c564b",
}

x = -df["Harmony#UMAP_Dimension_1"]
y = df["Harmony#UMAP_Dimension_2"]
c = df["majorclass"].map(pairwise_dict)
size = 0.5
# Create scatter plot
save = (
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/ATAC_majorclass_UMAP.tiff"
)

plot_umap(x, y, c, size, save)
