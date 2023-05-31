import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/AC_subtype_NRPC_obs.csv"
)
meta.index = meta["Unnamed: 0"].values
adata = adata[meta["Unnamed: 0"]]
adata.obs = meta

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
adata.var.highly_variable.loc[["PAX6", "TFAP2B", "GAD1", "GAD2", "SLC6A9"]] = True
adata.var.highly_variable.loc[["PAX6", "TFAP2B", "GAD1", "GAD2", "SLC6A9"]]
adata = adata[:, adata.var.highly_variable]
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pl.dotplot(
    adata,
    ["PAX6", "TFAP2B", "GAD1", "GAD2", "SLC6A9"],
    "subclass",
    dendrogram=True,
    swap_axes=True,
)
fig = plt.gcf()
fig.set_size_inches(10, 2)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/plot_AC_subgroup_DEG_dotplot.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
plt.clf()
