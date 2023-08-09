import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC.h5ad"
)

mv.velocity_embedding_stream(adata_result, basis="umap")
mv.velocity_graph(adata_result)
mv.latent_time(adata_result)
adata_result.obs.loc[:,"root_cells"] = 0
adata_result.obs.loc[(adata_result.obs.Time=="10w")&(adata_result.obs.Region=="Peripheral"),"root_cells"] = 1

mv.latent_time(adata_result)

adata_result.obs.to_csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_latent_time.csv")