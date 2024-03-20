import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/to_terminal_states.csv"
)
df.index = df["Unnamed: 0"].values
adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/NRPC.h5ad"
)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata, resolution=2)
adata = adata[adata.obs.index.isin(df.index)]

adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/NRPC_res2.h5ad"
)
