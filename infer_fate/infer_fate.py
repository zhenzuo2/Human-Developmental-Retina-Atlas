import scvelo as scv
import joblib
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
import numpy as np
import pandas as pd

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class_ldata_dynamic.h5ad",
    cache=False,
)

vk = joblib.load(
    "/storage/singlecell/zz4/fetal_bash/results/vk/merged_raw_filtered_umap_10000_major_sub_class_ldata_dynamic_vk.h5ad"
)
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.9 * vk + 0.1 * ck
g = GPCCA(combined_kernel)
print(g)

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "AC Precursor": "AC",
        "BC Precursor": "BC",
        "Cone Precursor": "Cone",
        "GABAergic": "AC",
        "Glycinergic": "AC",
        "HC0": "HC",
        "HC1": "HC",
        "MG": "MG",
        "ML_Cone": "Cone",
        "NRPC": "RPC",
        "OFF-BC": "BC",
        "OFF_MGC": "RGC",
        "ON-BC": "BC",
        "ON_MGC": "RGC",
        "PRPC": "RPC",
        "RBC": "BC",
        "RGC Precursor": "RGC",
        "Rod": "Rod",
        "Rod Precursor": "Rod",
        "SACs": "AC",
        "S_Cone": "Cone",
        "dual ACs": "AC",
    }
)

terminal_states = adata.obs["scpred_prediction"].replace("RPC", np.nan)
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities()

temp = pd.DataFrame(g.adata.obsm["to_terminal_states"])
temp.columns = list(g.adata.obsm["to_terminal_states"].names)
temp.index = g.adata.obs.index

temp.to_csv("/storage/singlecell/zz4/fetal_bash/results/vk/to_terminal_states.csv")

g.adata.write(
    "/storage/singlecell/zz4/fetal_bash/results/vk/merged_raw_filtered_umap_10000_major_sub_class_ldata_dynamic_g_adata.h5ad"
)