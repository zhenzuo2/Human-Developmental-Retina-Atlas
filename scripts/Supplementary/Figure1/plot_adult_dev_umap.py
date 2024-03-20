import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({"font.size": 12})

adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_wadult_umap_10000.h5ad"
)
category_mapping = {"Query": "Development", "Reference": "Adult"}

adata.obs["temp"] = adata.obs["batch"].map(category_mapping)
sc.pl.umap(
    adata,
    color="temp",
    size=0.5,
    legend_loc="right margin",
    title="",
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks([])
plt.yticks([])
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/overall_umap_with_adult_by_batch.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
matplotlib.rcParams.update({"font.size": 8})
sc.pl.umap(
    adata, color="leiden", size=0.5, legend_loc="on data", title="", frameon=False
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks([])
plt.yticks([])
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/overall_umap_with_adult_by_leiden.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
matplotlib.rcParams.update({"font.size": 12})
sc.pl.umap(
    adata,
    color="majorclass",
    size=0.5,
    legend_loc="on data",
    title="",
    palette={
        "MG": "#9467bd",
        "Rod": "#17becf",
        "BC": "#bcbd22",
        "RGC": "#d62728",
        "RPE": "#ff7f0e",
        "Cone": "#e377c2",
        "HC": "#2ca02c",
        "Microglia": "#1f77b4",
        "AC": "#8c564b",
    },
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks([])
plt.yticks([])
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/overall_umap_with_adult_by_majorclass.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)

adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
mapping_dict = {
    "PRPC":"MG",
    "NRPC":"MG",
}
adata.obs['majorclass'] = adata.obs['majorclass'].map(mapping_dict)
sc.pl.umap(
    adata,
    color="majorclass",
    size=0.5,
    legend_loc="on data",
    title="",
    palette={
        "MG": "#9467bd",
        "Rod": "#17becf",
        "BC": "#bcbd22",
        "RGC": "#d62728",
        "RPE": "#ff7f0e",
        "Cone": "#e377c2",
        "HC": "#2ca02c",
        "Microglia": "#1f77b4",
        "AC": "#8c564b",
    },
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks([])
plt.yticks([])
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/overall_umap_without_adult_by_majorclass.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
