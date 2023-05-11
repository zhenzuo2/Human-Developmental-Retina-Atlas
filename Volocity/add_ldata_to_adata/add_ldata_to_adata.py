import sys
import scvelo as scv
import anndata
import os

adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
ldata_file = "/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/combined.h5ad"
output_file = "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_ldata.h5ad"

try:
   os.makedirs(os.path.dirname(output_file))
except FileExistsError:
   # directory already exists
   pass

adata = scv.read(adata_file)
ldata = scv.read(ldata_file)

adata = scv.utils.merge(adata, ldata)
adata.write(output_file)