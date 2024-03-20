import os
import scvelo as scv
import joblib
import cellrank as cr
import numpy as np
import pandas as pd
import scanpy as sc

adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad",
    cache=False,
)
cells = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/AC_w_NRPC.csv")
adata = adata[cells["Unnamed: 0"].values]
scv.pp.filter_and_normalize(
    adata, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False
)

sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)
print(g)

terminal_states = adata.obs["subclass"].replace("NRPC", np.nan)
terminal_states = adata.obs["subclass"].replace('AC Precursor', np.nan)
g.set_terminal_states(terminal_states)
g.compute_fate_probabilities(n_jobs=1, backend='threading')

temp = pd.DataFrame(g.fate_probabilities)
temp.columns = list(g.fate_probabilities.names)
temp.index = g.adata.obs.index

temp.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/AC_infer_fate/AC_to_terminal_states.csv"
)

joblib.dump(
    g,
    "/storage/chentemp/zz4/adult_dev_compare/results/AC_infer_fate/AC_dynamics_g.pkl",
)