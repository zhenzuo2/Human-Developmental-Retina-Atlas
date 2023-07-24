import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

matplotlib.rcParams.update({"font.size": 25})
adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC.h5ad")
adata.raw = None
adata_result = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad")
adata[adata.obs.leiden.isin(["10"])].obs.to_csv("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_10.csv")

adata_result.obs["temp"] = ""
adata_result.obs["temp2"] = "MG fate"
adata_result.obs.loc[
    adata[adata.obs.leiden.isin(["10"])].obs.index, "temp"
] = "Neurogenesis fate"

time = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_10_monocle3_pseudotime.csv",sep = " ")
time["bin"] = pd.cut(time.x,2,labels = ["Early","Late"])

adata_result.obs.loc[
    time.loc[time.bin=="Early"].index, "temp2"
] = "Early Neurogenesis fate"
adata_result.obs.loc[
    time.loc[time.bin=="Late"].index, "temp2"
] = "Late Neurogenesis fate"
adata_result.obs.loc[adata[adata.obs.leiden.isin(["0"])].obs.index, "temp"] = "MG fate"
adata_result = adata_result[adata_result.obs.temp.isin(["Neurogenesis fate", "MG fate"])]
df = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/data/Human_TF_Full_List/TFsDatabaseExtract_v_1.01.csv")
TF = df.loc[df["Is TF?"]=="Yes","HGNC symbol"]
TF = [string for string in TF if not string.startswith("ZNF")]
sc.pp.filter_genes(adata_result, min_counts=1000)
sc.pp.normalize_total(adata_result)
sc.pp.log1p(adata_result)
adata_result = adata_result[:,adata_result.var.index.isin(TF)]
sc.tl.rank_genes_groups(adata_result, "temp")
markers0 = list(sc.get.rank_genes_groups_df(adata_result, group="MG fate").names[:10])
markers1 = list(sc.get.rank_genes_groups_df(adata_result, group="Neurogenesis fate").names[:10])

df = sc.get.rank_genes_groups_df(adata_result, group="Neurogenesis fate")
NRPC_genes = df.loc[df.logfoldchanges>0,"names"]

NRPC = adata_result[adata_result.obs.temp.isin(["Neurogenesis fate"])]
NRPC.obs["temp"] = "Late"
NRPC.obs.loc[
    time.loc[time.bin=="Early"].index, "temp"
] = "Early"
sc.tl.rank_genes_groups(NRPC, "temp")
df = sc.get.rank_genes_groups_df(NRPC, group="Early")
genes = df.loc[df.logfoldchanges>0,"names"]

markers2 = list([x for x in genes.values if x in NRPC_genes.values][:10])
df = sc.get.rank_genes_groups_df(NRPC, group="Late")
genes = df.loc[df.logfoldchanges>0,"names"]
markers3 = list([x for x in genes.values if x in NRPC_genes.values][:10])

plt.clf()
sc.pl.dotplot(
    adata_result,
    var_names={
        "MG fate enriched TFs": markers0,
        "Neurogenesis fate enriched TFs": markers1,
        "Early Neurogenesis fate enriched TFs": markers2,
        "Late Neurogenesis fate enriched TFs": markers3,
    },
    groupby="temp2",
    dendrogram=False,
)
plt.xticks(rotation=90)
plt.title("")
plt.xlabel("")
plt.ylabel("")
fig = plt.gcf()
fig.set_size_inches(15, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/RPC_fate.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
    backend="cairo",
)
plt.clf()

import seaborn as sns
adata_result.obs["WEEKS"] = adata_result.obs.Days/7
ax = sns.boxplot(data=adata_result.obs, x="temp", y="WEEKS",hue = "Region")
ax.set_xlabel("")
ax.set_ylabel("PCW")
plt.setp(ax.get_xticklabels(), rotation=45)
plt.legend(bbox_to_anchor=(1.02, 0.1), loc='upper left', borderaxespad=0)
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/RPC_fate_clock_time.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
    backend="cairo",
)
plt.clf()


adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_MG_NRPC.h5ad")
NRPC = adata[NRPC.obs.index]
NRPC.obs["time"]=time.loc[NRPC.obs.index,"x"]
sc.pl.umap(NRPC,color = "time",frameon=False,size =20)
plt.title("")
plt.xlabel("")
plt.ylabel("")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/NRPC_fate_pseudo_time.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
    backend="cairo",
)
plt.clf()