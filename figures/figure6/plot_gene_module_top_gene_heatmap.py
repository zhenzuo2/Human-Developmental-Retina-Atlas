import joblib
import scvelo as scv
import plotly.express as px
import scanpy as sc
import hotspot
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib
from scipy.cluster.hierarchy import leaves_list
matplotlib.rcParams.update({'font.size': 12})

output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"
hs = joblib.load(
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/latent_time/PPRC_hs.create_modules.pkl"
)
adata = joblib.load(
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/latent_time/PPRC_hs.adata.pkl"
)

modules = hs.create_modules(min_gene_threshold=120, core_only=True, fdr_threshold=0.05)
hs.modules = hs.modules.replace({4: 1, 3: 2, 2: 3, 5: 4, 1: 5})

ii = leaves_list(hs.linkage)

mod_reordered = hs.modules.iloc[ii]

mod_map = {}
y = np.arange(modules.size)

for x in mod_reordered.unique():
    if x == -1:
        continue

mod_map[x] = y[mod_reordered == x].mean()

mod_reordered.to_csv("/storage/singlecell/zz4/fetal_snakemake/figures/figure6/mod_reordered.csv")
adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
PRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC/adata_umap.h5ad"
)
adata_result = adata_result[PRPC.obs.index]
adata_result.obsm["X_umap"] = PRPC[adata_result.obs.index].obsm["X_umap"]

df = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_latent_time.csv"
)
df.index = df["Unnamed: 0"].values
adata_result.obs["latent_time"]=df.loc[adata_result.obs.index,"latent_time"].values

sc.pp.normalize_total(adata_result)
sc.pp.scale(adata_result)

df = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/data/Human_TF_Full_List/TFsDatabaseExtract_v_1.01.csv")
TF = df.loc[df["Is TF?"]=="Yes","HGNC symbol"]

mod_reordered= mod_reordered[mod_reordered.index.isin(TF)]
top_genes = mod_reordered[mod_reordered==5][:10].index
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='latent_time',yticklabels=True,vmax= 0.5)
fig = plt.gcf()
fig.set_size_inches(10,2)
plt.savefig(
    output_file_path + "PRPC_gene_module_5_heatmap.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

top_genes = mod_reordered[mod_reordered==4][:10].index
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='latent_time',yticklabels=True,vmax= 0.5)
fig = plt.gcf()
fig.set_size_inches(10,2)
plt.savefig(
    output_file_path + "PRPC_gene_module_4_heatmap.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

top_genes = mod_reordered[mod_reordered==3][:10].index
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='latent_time',yticklabels=True,vmax= 0.6)
fig = plt.gcf()
fig.set_size_inches(10,2)
plt.savefig(
    output_file_path + "PRPC_gene_module_3_heatmap.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

top_genes = mod_reordered[mod_reordered==2][:10].index
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='latent_time',yticklabels=True,vmax= 0.5)
fig = plt.gcf()
fig.set_size_inches(10,2)
plt.yticks(rotation = 0)
plt.savefig(
    output_file_path + "PRPC_gene_module_2_heatmap.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

top_genes = mod_reordered[mod_reordered==1][:10].index
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='latent_time',yticklabels=True,vmax= 0.2)
fig = plt.gcf()
fig.set_size_inches(10,2)
plt.savefig(
    output_file_path + "PRPC_gene_module_1_heatmap.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)