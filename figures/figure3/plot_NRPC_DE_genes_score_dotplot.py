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
NRPC.obs["subclass"] = "Undetermined"

sc.pp.neighbors(NRPC, use_rep="X_scVI")
sc.tl.leiden(NRPC, resolution=3)

clusters = ["28", "3", "13", "15"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["6", "30","35","25"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["21","14","27"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["33","0"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["12","22"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["4", "19","23"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"
NRPC.obs.to_csv("/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/Annotated_NRPC.csv")
sc.pp.highly_variable_genes(NRPC,flavor='seurat_v3',n_top_genes=5000,subset=True)
sc.pp.normalize_total(NRPC)
sc.pp.log1p(NRPC)

gs = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/gene_score.h5ad",
)
sc.pp.normalize_total(gs)
gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]

cells = [x for x in NRPC.obs.index if x in gs.obs.index]
NRPC = NRPC[cells]
gs = gs[cells]

gs.obs['subclass'] = NRPC.obs['subclass']
gs = gs[:,[s for s in gs.var.index if not s.startswith("MIR")]]

gs.obs["subclass"] = pd.Categorical(
    list(gs.obs["subclass"]),
    categories=["RGC","AC","HC","Rod","Undetermined","BC","Cone"],
)

sc.tl.rank_genes_groups(gs, "subclass")
sc.pl.rank_genes_groups_dotplot(gs, n_genes=10,cmap="YlGnBu",groups = ["RGC","AC","HC","Rod","Undetermined","BC","Cone"])
fig = plt.gcf()
fig.set_size_inches(20, 4)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_DE_genes_score_dotplot.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()