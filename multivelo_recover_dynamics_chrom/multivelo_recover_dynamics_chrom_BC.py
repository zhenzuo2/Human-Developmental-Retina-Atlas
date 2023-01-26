import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys

adata_rna_file = "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_run_umap_BC/adata_umap.h5ad"
adata_atac_file = (
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_knn_smooth_chrom/BC.h5ad"
)
output_file = "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/BC.h5ad"

adata_rna = scv.read(adata_rna_file)
All = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_modules_meta_All_feature_selection_FALSE.csv"
)
Macula = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_modules_meta_Macula_feature_selection_FALSE.csv"
)
Peripheral = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_modules_meta_Peripheral_feature_selection_FALSE.csv"
)
target = list(
    set(list(set(All.target)) + list(set(Macula.target)) + list(set(Peripheral.target)))
)
tf = list(set(list(set(All.tf)) + list(set(Macula.tf)) + list(set(Peripheral.tf))))
tf = list(set(tf))
adata_rna = adata_rna[:, [x for x in adata_rna.var.index if x in target + tf]]

adata_atac = scv.read(adata_atac_file)

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
    n_jobs=20,
)
adata_result.write(output_file)