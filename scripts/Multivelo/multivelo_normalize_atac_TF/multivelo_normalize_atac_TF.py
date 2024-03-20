import scvelo as scv
import multivelo as mv
import scanpy as sc
import sys
import os

input_file_path = sys.argv[1]
sample = sys.argv[2]
output_result_path = sys.argv[3]

try:
   os.makedirs(output_result_path)
except FileExistsError:
   # directory already exists
   pass

input_path = sys.argv[1]

adata_atac = scv.read(input_file_path)

sc.pp.filter_cells(adata_atac, min_counts=1000)
sc.pp.filter_cells(adata_atac, max_counts=40000)

# We normalize aggregated peaks with TF-IDF.
mv.tfidf_norm(adata_atac)

adata_atac.write(output_result_path + sample + ".h5ad")