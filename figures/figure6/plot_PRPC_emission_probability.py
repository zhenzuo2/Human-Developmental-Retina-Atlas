import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
from cellrank.tl.kernels import VelocityKernel
import joblib
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
import seaborn as sns

output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_MG_NRPC.h5ad"
)

scv.tl.velocity(adata_result)
vk = VelocityKernel(adata_result)
vk.compute_transition_matrix()

ck = ConnectivityKernel(adata_result).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = GPCCA(combined_kernel)

print(g)

terminal_states = adata_result.obs["majorclass"].replace("PRPC", np.nan)
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities()

temp = pd.DataFrame(g.adata.obsm["to_terminal_states"])
temp.columns = list(g.adata.obsm["to_terminal_states"].names)
temp.index = g.adata.obs.index
adata_result.obs["MG"] = temp.MG
adata_result.obs["NRPC"] = temp.NRPC

adata_result.write("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_MG_NRPC_TO_MG_NRPC_Fate.h5ad")

sc.pl.umap(adata_result[adata_result.obs.majorclass=="PRPC",],color = "NRPC",vmax = 0.86,vmin = 0.7,frameon=False,title = "",colorbar_loc = None,size=50)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_emission_p_to_NRPC.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

sc.pl.umap(adata_result[adata_result.obs.majorclass=="PRPC",],color = "MG",vmax = 0.8,vmin = 0.6,frameon=False,title = "",colorbar_loc = None,size=50)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_emission_p_to_MG.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)



df = adata_result.obs.loc[adata_result.obs.leiden.isin(['0',
 '1',
 '10',
 #'11',
 #'12',
 #'13',
 '14',
 '15',
 '16',
 '17',
 '2',
 '3',
 '4',
 #'5',
 #'6',
 '7',
 '8',
 '9']),:]

sns.boxplot(data=df, x="leiden", y="MG")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_emission_p_to_MG_boxplot.svg",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)
plt.clf()

sns.boxplot(data=df, x="leiden", y="NRPC")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_emission_p_to_NRPC_boxplot.svg",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)
plt.clf()

sc.pl.umap(adata_result,color = "leiden",legend_loc='on data',size=20,frameon = False)
fig = plt.gcf()
fig.set_size_inches(6, 6)
plt.savefig(
    output_file_path + "PRPC_leiden.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)
plt.clf()