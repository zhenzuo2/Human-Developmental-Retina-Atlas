import os
import scvelo as scv
import joblib
import cellrank as cr
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

g = joblib.load(
    "/storage/chentemp/zz4/adult_dev_compare/results/PRPC_infer_fate/PRPC_dynamics_g.pkl"
)

adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad",
    cache=False,
)
adata = adata[adata.obs.majorclass.isin(["PRPC"])]
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

driver_clusters = ["MG", "NRPC"]
delta_df = g.compute_lineage_drivers(
    lineages=["NRPC"], cluster_key="majorclass", clusters=driver_clusters
)
g.adata.obs["Fate Probabilities NRPC"] = g.fate_probabilities["NRPC"].X.flatten()

adata.obs["Fate Probabilities NRPC"] = g.adata[adata.obs.index].obs[
    "Fate Probabilities NRPC"
]
df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_Module_gene_score.csv")
df.index = df['Unnamed: 0'].values
adata = adata[adata.obs.index.isin(df.index)]
adata = adata[df.index,:]
adata.obs["Module1"] = df.loc[adata.obs.index,"Module1"]
adata.obs["Module2"] = df.loc[adata.obs.index,"Module2"]
adata.obs["Module3"] = df.loc[adata.obs.index,"Module3"]

from scipy.stats import pearsonr

def calculate_pearson_correlation(x, y):
    """
    Calculate Pearson correlation coefficient and p-value between two arrays.

    Parameters:
    - x, y: Input arrays for which to calculate the correlation.

    Returns:
    - pearson_corr: Pearson correlation coefficient.
    - p_value: Two-tailed p-value.
    """
    pearson_corr, p_value = pearsonr(x, y)
    return pearson_corr, p_value

corr_coefficient, p_value = calculate_pearson_correlation(adata.obs["Fate Probabilities NRPC"], adata.obs["Module1"])
print(f"Pearson Correlation Coefficient: {corr_coefficient}")
print(f"P-value: {p_value}")



