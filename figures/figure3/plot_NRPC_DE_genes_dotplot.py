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

clusters = ["2", "13"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["0", "1", "18"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["9", "10"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["4"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["17"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["7", "14"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"

sc.pp.highly_variable_genes(NRPC,flavor='seurat_v3',n_top_genes=5000,subset=True)
sc.pp.normalize_total(NRPC)
sc.pp.log1p(NRPC)

sc.tl.rank_genes_groups(NRPC, "subclass")
sc.pl.rank_genes_groups_dotplot(NRPC, n_genes=10)
fig = plt.gcf()
fig.set_size_inches(20, 4)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_DE_genes_dotplot.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()