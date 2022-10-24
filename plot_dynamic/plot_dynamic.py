import scvelo as scv
import scanpy as sc
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

scv.settings.plot_prefix = ""

adata = scv.read(input_path)

scv.pl.velocity_embedding_stream(adata, basis='umap', save = output_path)