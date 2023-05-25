# Import packages
# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns

sc.set_figure_params(transparent=True, fontsize=10)
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

gs = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/gene_score.h5ad",
)

gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]

cells = [x for x in NRPC.obs.index if x in gs.obs.index]
NRPC = NRPC[cells]
gs = gs[cells]

gs.obs['subclass'] = NRPC.obs['subclass']
gs = gs[:,[s for s in gs.var.index if not s.startswith("MIR")]]

gs.obs["subclass"] = pd.Categorical(
    list(gs.obs["subclass"]),
    categories=["RGC","AC","HC","Rod","Unknown","BC","Cone"],
)

sc.tl.rank_genes_groups(gs, "subclass")
sc.pl.rank_genes_groups_dotplot(gs, n_genes=10,cmap="YlGnBu",groups = ["RGC","AC","HC","Rod","Unknown","BC","Cone"])
fig = plt.gcf()
fig.set_size_inches(20, 4)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_DE_genes_score_dotplot.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()