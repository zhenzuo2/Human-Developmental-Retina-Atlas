import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad"
)
adata = adata[adata.obs.majorclass == "PRPC"]

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata, n_jobs=10)

scv.tl.velocity_pseudotime(adata)

adata.obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/infer_velocity_pseudotime/PRPC.csv"
)
