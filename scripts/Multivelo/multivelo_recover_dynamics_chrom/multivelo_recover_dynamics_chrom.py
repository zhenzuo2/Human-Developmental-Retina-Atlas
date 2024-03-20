import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

samples = sys.argv[1]
pd.options.display.max_columns = None
scv.set_figure_params(dpi=600, dpi_save=600)
input_dir = "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_run_umap"
adata_rna_file = input_dir + "/" + samples + ".h5ad"
adata_atac_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_knn_smooth_chrom/"
    + samples
    + ".h5ad"
)
output_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/"
    + samples
    + ".h5ad"
)

adata_rna = scv.read(adata_rna_file)

sc.pp.highly_variable_genes(
    adata_rna, flavor="seurat_v3", n_top_genes=2000, subset=True
)
adata_atac = scv.read(adata_atac_file)

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna)

print(adata_rna)
print(adata_atac)

sc.pl.umap(adata_rna, color="majorclass")

adata_result = mv.recover_dynamics_chrom(
    adata_rna,
    adata_atac,
    max_iter=5,
    init_mode="invert",
    verbose=True,
    parallel=True,
    save_plot=False,
    rna_only=False,
    fit=True,
    n_anchors=500,
    extra_color_key="majorclass",
    n_jobs=10,
)
adata_result.write(output_file)