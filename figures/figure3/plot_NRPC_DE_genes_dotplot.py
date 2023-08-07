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
clusters = ["33","0"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["12","22"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["4", "19","23"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"

sc.pp.highly_variable_genes(NRPC,flavor='seurat_v3',n_top_genes=5000,subset=True)
sc.pp.normalize_total(NRPC)
sc.pp.log1p(NRPC)

df = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/data/Human_TF_Full_List/TFsDatabaseExtract_v_1.01.csv")
TF = df.loc[df["Is TF?"]=="Yes","HGNC symbol"]
NRPC= NRPC[:,NRPC.var.index.isin(TF)]
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


adata = NRPC[NRPC.obs.leiden.isin(['27','14',"21"])]
adata.obs['leiden'] = adata.obs.leiden.map({'21': 'TBD', '27': 'OFF Cone BC','14':"ON Cone BC/RBC"})
sc.tl.rank_genes_groups(adata, "leiden")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=10)
fig = plt.gcf()
fig.set_size_inches(20, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/BC_subtype_DE_genes_dotplot.svg",
    dpi=600,
    bbox_inches="tight",
)
