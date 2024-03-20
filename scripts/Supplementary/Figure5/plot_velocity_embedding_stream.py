import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt


directory_path = "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/"

# Get the list of files in the directory
files_list = os.listdir(directory_path)

for f in files_list:
    adata_result = sc.read_h5ad(directory_path + f)
    label = os.path.splitext(f)[0]
    scv.pp.remove_duplicate_cells(adata_result)

    mv.velocity_embedding_stream(adata_result)
    adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
    mv.velocity_embedding_stream(
        adata_result,
        color="majorclass",
        title="",
        alpha=1,
        palette={
            "MG": "#9467bd",
            "Rod": "#17becf",
            "BC": "#bcbd22",
            "RGC": "#d62728",
            "NRPC": "#ff7f0e",
            "Cone": "#e377c2",
            "HC": "#2ca02c",
            "PRPC": "#1f77b4",
            "AC": "#8c564b",
        },
    )
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure5/"
        + label
        + "_velocity_embedding_stream.tiff",
        dpi=300,
        transparent=True,
    )
    mv.velocity_embedding_stream(adata_result, color="Days", title="", alpha=1)
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure5/"
        + label
        + "_velocity_embedding_stream_Days.tiff",
        dpi=300,
        transparent=True,
    )
#  For Cone

adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
f = "Cone.h5ad"
adata_result = sc.read_h5ad(directory_path + f)
label = os.path.splitext(f)[0]
scv.pp.remove_duplicate_cells(adata_result)

mv.velocity_embedding_stream(adata_result)
adata_result.obsm["X_umap"] = adata[adata_result.obs.index].obsm["X_umap"]
adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
sc.tl.leiden(adata_result)
sc.pl.umap(adata_result, color="leiden", legend_loc="on data")
#adata_result = adata_result[adata_result.obs.leiden != "1"]
mv.velocity_embedding_stream(
    adata_result,
    color="majorclass",
    title="",
    alpha=1,
    palette={
        "MG": "#9467bd",
        "Rod": "#17becf",
        "BC": "#bcbd22",
        "RGC": "#d62728",
        "NRPC": "#ff7f0e",
        "Cone": "#e377c2",
        "HC": "#2ca02c",
        "PRPC": "#1f77b4",
        "AC": "#8c564b",
    },
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure5/"
    + label
    + "_velocity_embedding_stream.tiff",
    dpi=300,
    transparent=True,
)

mv.velocity_embedding_stream(adata_result, color="Days", title="", alpha=1)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure5/"
    + label
    + "_velocity_embedding_stream_Days.tiff",
    dpi=300,
    transparent=True,
)



adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
f = "Rod.h5ad"
adata_result = sc.read_h5ad(directory_path + f)
label = os.path.splitext(f)[0]
scv.pp.remove_duplicate_cells(adata_result)

mv.velocity_embedding_stream(adata_result)
adata_result.obsm["X_umap"] = adata[adata_result.obs.index].obsm["X_umap"]
adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
sc.tl.leiden(adata_result)
sc.pl.umap(adata_result, color="leiden", legend_loc="on data")
mv.velocity_embedding_stream(
    adata_result,
    color="majorclass",
    title="",
    alpha=1,
    palette={
        "MG": "#9467bd",
        "Rod": "#17becf",
        "BC": "#bcbd22",
        "RGC": "#d62728",
        "NRPC": "#ff7f0e",
        "Cone": "#e377c2",
        "HC": "#2ca02c",
        "PRPC": "#1f77b4",
        "AC": "#8c564b",
    },
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure5/"
    + label
    + "_velocity_embedding_stream.tiff",
    dpi=300,
    transparent=True,
)

mv.velocity_embedding_stream(adata_result, color="Days", title="", alpha=1)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure5/"
    + label
    + "_velocity_embedding_stream_Days.tiff",
    dpi=300,
    transparent=True,
)
