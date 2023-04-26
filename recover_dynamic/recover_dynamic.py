import scvelo as scv
import scanpy as sc
import sys
import os

input_path = sys.argv[1]
output_path = sys.argv[2]

try:
   os.makedirs(os.path.dirname(output_path))
except FileExistsError:
   # directory already exists
   pass

def velocity_graph_stochastic(input_path):
    print("Runing stochastic for " + input_path)
    adata = scv.read(input_path)
    adata = adata[adata.obs.Region.isin(['Macula','Peripheral'])&(adata.obs['Data Type']=='Multiomics')]
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)
    
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata, n_jobs=20)
    return adata


velocity_graph_stochastic(
    input_path
).write(
    output_path
)