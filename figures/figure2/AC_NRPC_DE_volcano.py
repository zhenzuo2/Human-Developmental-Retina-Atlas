import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/AC.h5ad")
adata_result = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad")
adata_result = adata_result[adata.obs.index[adata.obs.leiden=="9"]]
sc.pp.highly_variable_genes(adata_result,flavor='seurat_v3',n_top_genes=2000,subset=True)
sc.pp.normalize_total(adata_result)
sc.pp.log1p(adata_result)
sc.pp.neighbors(adata_result)
sc.tl.umap(adata_result)
sc.pl.umap(adata_result,color = "Days")
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_Days.svg",
    dpi=600,
    bbox_inches="tight",
)
sc.tl.leiden(adata_result,resolution=0.5)
sc.pl.umap(adata_result,color = "leiden",legend_loc = "on data")
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_leiden.svg",
    dpi=600,
    bbox_inches="tight",
)
sc.tl.rank_genes_groups(adata_result, "leiden")
sc.get.rank_genes_groups_df(adata_result,group="4").to_csv("/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_DE.csv")