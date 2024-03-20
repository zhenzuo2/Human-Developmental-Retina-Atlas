import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/NRPC_res2.h5ad"
)
adata.obs["subclass"] = adata.obs["subclass"].astype(str)
clusters = ["0"]
adata.obs.loc[adata.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["15"]
adata.obs.loc[adata.obs.leiden.isin(clusters), "subclass"] = "RGC"
clusters = ["5"]
adata.obs.loc[adata.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["1"]
adata.obs.loc[adata.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["25","2"]
adata.obs.loc[adata.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["4"]
adata.obs.loc[adata.obs.leiden.isin(clusters), "subclass"] = "BC"

adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.h5ad"
)
adata.obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.csv"
)
