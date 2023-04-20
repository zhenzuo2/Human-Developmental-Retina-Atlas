import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys

adata_rna_file = sys.argv[1]
adata_atac_file = sys.argv[2]
output_file = sys.argv[3]

adata_rna = scv.read(adata_rna_file)
adata_atac = scv.read(adata_atac_file)

sc.pp.highly_variable_genes(
        adata_rna, flavor="seurat_v3", n_top_genes=1000, subset=True
    )

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)

print(adata_rna)
print(adata_atac)

sc.pl.umap(adata_rna, color="majorclass")

adata_result = mv.recover_dynamics_chrom(
    adata_rna,
    adata_atac,
    max_iter=5,
    init_mode="invert",
    verbose=False,
    parallel=True,
    save_plot=False,
    rna_only=False,
    fit=True,
    n_anchors=500,
    extra_color_key="majorclass",
    n_jobs=8,
)

# Save the result for use later on
adata_result.write(output_file)