import scanpy as sc
import scmer

adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

sc.pp.subsample(adata, n_obs=10000)

from scmer import UmapL1
model_20 = UmapL1.tune(target_n_features=300, X=adata.X.toarray(),n_threads=10,max_iter=2,max_outer_iter=2)
print(*adata.var_names[model_20.get_mask()])




