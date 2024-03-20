import scanpy as sc

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_expression.h5ad"
)