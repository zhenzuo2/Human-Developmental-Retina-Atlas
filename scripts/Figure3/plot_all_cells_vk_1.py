import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import joblib

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2
adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad"
)

scv.pp.filter_and_normalize(
    adata, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False
)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="stochastic")

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
vk.plot_projection(color="majorclass")
joblib.dump(
    vk,
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata_dynamics_vk.pkl"
)

