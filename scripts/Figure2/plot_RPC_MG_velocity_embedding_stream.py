import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt


adata_result = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/RPC_MG.h5ad"
)
mv.velocity_embedding_stream(adata_result)

adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
mv.velocity_embedding_stream(adata_result, color="Days",title = "",alpha = 1)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/RPC_MG_mv.velocity_embedding_stream.Days.png",
    dpi=300,
    transparent=True,
)
