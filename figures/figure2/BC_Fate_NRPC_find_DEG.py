# Import packages
# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns

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
clusters = ["14"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cluster 2"
clusters = ["27"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cluster 2"
clusters = ["21"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cluster 1"
clusters = ["33","0"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["12","22"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cluster 3"
clusters = ["4", "19","23"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"


adata = NRPC[NRPC.obs.subclass.isin(["Cluster 1","Cluster 3"])]
sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", n_top_genes=5000, subset=True
)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, 'subclass')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=15)
fig = plt.gcf()
fig.set_size_inches(15, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/BC_Fate_NRPC_find_DEG.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()