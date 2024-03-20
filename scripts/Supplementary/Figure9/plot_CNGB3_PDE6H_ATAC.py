import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_score_imputed_umap.h5ad"
)

sc.pl.umap(
    adata[adata.obs.Region == "Macula"], color="CNGB3", s=10, frameon=False, vmax=0.5
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_CNGB3_M_gs.tiff",
    transparent=True,
    dpi=300,
)

sc.pl.umap(
    adata[adata.obs.Region == "Peripheral"],
    color="CNGB3",
    s=10,
    frameon=False,
    vmax=0.5
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_CNGB3_P_gs.tiff",
    transparent=True,
    dpi=300,
)



sc.pl.umap(
    adata[adata.obs.Region == "Macula"], color="PDE6H", s=10, frameon=False, vmax=0.5
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_PDE6H_M_gs.tiff",
    transparent=True,
    dpi=300,
)

sc.pl.umap(
    adata[adata.obs.Region == "Peripheral"],
    color="PDE6H",
    s=10,
    frameon=False,
    vmax=0.5
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure9/RNA_PDE6H_P_gs.tiff",
    transparent=True,
    dpi=300,
)