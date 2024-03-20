import scanpy as sc

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

adata.write("/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_expression.h5ad")