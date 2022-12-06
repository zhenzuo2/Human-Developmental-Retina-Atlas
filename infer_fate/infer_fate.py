import scvelo as scv
import joblib
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
import numpy as np
import pandas as pd

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/merged_h5ad_adult_annotated_umap.h5ad",
    cache=False,
)

vk = joblib.load(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/merged_h5ad_adult_annotated_umap_vk.pkl"
)
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.9 * vk + 0.1 * ck
g = GPCCA(combined_kernel)
print(g)

terminal_states = adata.obs["scpred_prediction"].replace("RPC", np.nan)
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities()

temp = pd.DataFrame(g.adata.obsm["to_terminal_states"])
temp.columns = list(g.adata.obsm["to_terminal_states"].names)
temp.index = g.adata.obs.index

temp.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/to_terminal_states.csv"
)

g.adata.write(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/merged_h5ad_adult_annotated_umap_g_adata.h5ad"
)