import scanpy as sc

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad")

adata.obs.to_csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/ALL.csv")