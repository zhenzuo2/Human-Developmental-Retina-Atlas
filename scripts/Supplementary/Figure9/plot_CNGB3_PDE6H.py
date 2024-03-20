import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

sc.pl.umap(
    adata[adata.obs.Region == "Macula"], color="CNGB3", s=10, frameon=False, vmax=3.5
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_CNGB3_M.tiff",
    transparent=True,
    dpi=300,
)

sc.pl.umap(
    adata[adata.obs.Region == "Peripheral"],
    color="CNGB3",
    s=10,
    frameon=False,
    vmax=3.5,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_CNGB3_P.tiff",
    transparent=True,
    dpi=300,
)


sc.pl.umap(
    adata[adata.obs.Region == "Macula"],
    color="PDE6H",
    s=10,
    vmax=3.5,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_PDE6H_M.tiff",
    transparent=True,
    dpi=300,
)

sc.pl.umap(
    adata[adata.obs.Region == "Peripheral"],
    color="PDE6H",
    s=10,
    vmax=3.5,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_PDE6H_P.tiff",
    transparent=True,
    dpi=300,
)

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_score_imputed_umap.h5ad"
)
sc.pl.umap(
    adata[adata.obs.Region == "Macula"], color="CNGB3", s=10, frameon=False, vmax=0.4
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/ATAC_CNGB3_M.tiff",
    transparent=True,
    dpi=300,
)

sc.pl.umap(
    adata[adata.obs.Region == "Peripheral"],
    color="CNGB3",
    s=10,
    frameon=False,
    vmax=0.4,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/ATAC_CNGB3_P.tiff",
    transparent=True,
    dpi=300,
)


sc.pl.umap(
    adata[adata.obs.Region == "Macula"], color="PDE6H", s=10, frameon=False, vmax=0.4
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/ATAC_PDE6H_M.tiff",
    transparent=True,
    dpi=300,
)

sc.pl.umap(
    adata[adata.obs.Region == "Peripheral"],
    color="PDE6H",
    s=10,
    frameon=False,
    vmax=0.4,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/ATAC_PDE6H_P.tiff",
    transparent=True,
    dpi=300,
)
