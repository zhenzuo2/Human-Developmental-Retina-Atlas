# Import packages
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read(
"/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
sc.pl.umap(
    adata,
    color="celltype",
    size=2,
    title="",
    frameon=False,
    legend_loc="on data",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/cell_type_annotation.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)