import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 20})

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata = adata[adata.obs.Region.isin(["Macula", "Peripheral"])]
adata = adata[adata.obs.majorclass == "Cone"]

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

sc.pl.violin(adata, keys=["CNGB3"], groupby="Region")
plt.xlabel("")
plt.ylabel("")
plt.title("CNGB3")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure7/FovealHypoplasia_violin_CNGB3.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.violin(adata, keys=["PDE6H"], groupby="Region")
plt.xlabel("")
plt.ylabel("")
plt.title("PDE6H")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure7/FovealHypoplasiaHetamap_violin_PDE6H.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
