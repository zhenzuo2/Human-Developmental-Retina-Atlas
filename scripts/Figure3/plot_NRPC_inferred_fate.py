import scanpy as sc
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["scatter.edgecolors"] = "none"
font = {"family": "Arial", "size": 12}
plt.rc("font", **font)

def plot_umap(x, y, c, size, save):
    plt.figure(figsize=(5, 5))
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


adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.h5ad"
)
sc.pl.umap(
    adata,
    size=20,
    legend_loc=None,
    colorbar_loc=None,
    color="leiden",
    frameon=False,
    title = ""
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_leiden.tiff",
    dpi=300,
    bbox_inches="tight",
)
plt.clf()

pairwise_dict = {
    "MG": "#9467bd",
    "Rod": "#17becf",
    "BC": "#bcbd22",
    "RGC": "#d62728",
    "NRPC": "#d3d3d3",
    "Cone": "#e377c2",
    "HC": "#2ca02c",
    "PRPC": "#1f77b4",
    "AC": "#8c564b",
}

x = [x[0] for x in adata.obsm["X_umap"]]
y = [x[1] for x in adata.obsm["X_umap"]]

c = adata.obs["subclass"].map(pairwise_dict)
size = 3
# Create scatter plot
save = (
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC_inferred_fate.tiff"
)
plot_umap(x, y, c, size, save)
