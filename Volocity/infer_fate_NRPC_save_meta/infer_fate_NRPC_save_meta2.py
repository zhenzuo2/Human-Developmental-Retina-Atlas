# Import packages
# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)
sc.set_figure_params(transparent=True, fontsize=15)
NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC = NRPC[NRPC.obs.leiden!="12"]
NRPC.obs["subclass"] = "Unknown"

sc.pp.neighbors(NRPC, use_rep="X_scVI")
sc.tl.leiden(NRPC, resolution=3)

clusters = ["28", "3", "13", "15"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["6", "30","35","25"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["21","14","27"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["33","0"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["12","22"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["4", "19","23"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"

clusters = ["21","14","27"]
pd.concat(
    [
        adata[NRPC[(NRPC.obs.leiden.isin(clusters))].obs.index].obs,
        adata[adata.obs.scpred_prediction == "BC"].obs,
    ]
).to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_BC.csv"
)