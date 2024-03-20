import scvelo as scv
import scanpy as sc
import sys
import os

input_path = "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad"
output_path = "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata_dynamics.h5ad"

try:
    os.makedirs(os.path.dirname(output_path))
except FileExistsError:
    # directory already exists
    pass


def scv_recover_dynamics(input_path):
    print("Runing stochastic for " + input_path)
    adata = scv.read(input_path)
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=10000)
    scv.pp.log1p(adata)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata)
    return adata


scv_recover_dynamics(input_path).write(output_path)
