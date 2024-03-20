import os
os.environ["NUMBA_DISABLE_JIT"] = "1"
import scvelo as scv
import joblib
import cellrank as cr
import numpy as np
import pandas as pd


adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata_dynamics.h5ad",
    cache=False,
)

vk = joblib.load(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata_dynamics_vk.h5ad"
)
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)
print(g)

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)

terminal_states = adata.obs["scpred_prediction"].replace("RPC", np.nan)
g.set_terminal_states(terminal_states)
g.compute_fate_probabilities(n_jobs=1, backend='threading')

temp = pd.DataFrame(g.fate_probabilities)
temp.columns = list(g.fate_probabilities.names)
temp.index = g.adata.obs.index

temp.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/to_terminal_states.csv"
)

joblib.dump(
    g,
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata_dynamics_g.h5ad",
)