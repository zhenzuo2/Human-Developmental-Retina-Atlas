import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

adata_result = sc.read_h5ad("/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/RPC_MG.h5ad")
mv.velocity_embedding_stream(adata_result)

adata_result.obs['Days'] = adata_result.obs['Days'].astype(float)

mv.velocity_graph(adata_result,n_jobs = 10)
mv.latent_time(adata_result)


adata_result.obs.to_csv("/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_latent_time/PRPC.csv")