import sys
import scvelo as scv
import anndata
import os

adata_file = "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
ldata_file = "/storage/chentemp/zz4/adult_dev_compare/results/RNA_velocity/combined.h5ad"
output_file = "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad"

try:
   os.makedirs(os.path.dirname(output_file))
except FileExistsError:
   # directory already exists
   pass

adata = scv.read(adata_file)
ldata = scv.read(ldata_file)

adata = scv.utils.merge(adata, ldata)
adata.write(output_file)