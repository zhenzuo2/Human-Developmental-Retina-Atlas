import scvelo as scv
import joblib
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
import numpy as np
import pandas as pd

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata_dynamics.h5ad",
    cache=False,
)

vk = joblib.load(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata_dynamics_vk.h5ad"
)
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = GPCCA(combined_kernel)
print(g)

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)

terminal_states = adata.obs["scpred_prediction"].replace("RPC", np.nan)
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities()

temp = pd.DataFrame(g.adata.obsm["to_terminal_states"])
temp.columns = list(g.adata.obsm["to_terminal_states"].names)
temp.index = g.adata.obs.index

temp.to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/to_terminal_states.csv"
)

joblib.dump(
    g,
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata_dynamics_g.h5ad",
)
