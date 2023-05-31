import scanpy as sc
import pandas as pd

adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/AC_subtype_NRPC_obs.csv"
)
meta.index = meta["Unnamed: 0"].values
adata = adata[meta["Unnamed: 0"]]
adata.obs = meta
adata = adata[adata.obs.majorclass == "Glycinergic"]
adata.obs["temp"] = adata.obs.subclass.isin(["AC4", "AC10"])
adata.obs["temp"] = adata.obs["temp"].astype(str)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000, subset=True)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, "temp")
sc.get.rank_genes_groups_df(adata, group="False").to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/AC_subtype_DEGs.csv"
)
