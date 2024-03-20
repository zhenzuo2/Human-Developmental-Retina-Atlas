import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

matplotlib.rcParams.update({"font.size": 25})
adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC_MG_NRPC/adata_umap.h5ad")
adata.obs["temp"] = "Rest of RPC"
adata = adata[adata.obs.majorclass.isin(["PRPC","NRPC"])]
adata.obs.loc[adata.obs.leiden.isin(["14"])&(adata.obs.majorclass=="PRPC"), "temp"] = "Neurogenesis"
#adata.obs.loc[adata.obs.leiden.isin(["0","17"])&(adata.obs.majorclass=="PRPC"), "temp"] = "MG precursor"
adata.obs.loc[adata.obs.majorclass=="NRPC", "temp"] = "NRPC"
adata[adata.obs.leiden.isin(["14"])&(adata.obs.majorclass=="PRPC")].obs.to_csv("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_10.csv")

adata_result = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad")
adata_result = adata_result[adata.obs.index,:]
adata_result.obs['temp'] = adata.obs['temp']

sc.pp.filter_genes(adata_result, min_counts=1000)
sc.pp.normalize_total(adata_result)
sc.pp.log1p(adata_result)
sc.tl.rank_genes_groups(adata_result, "temp")
#markers0 = list(sc.get.rank_genes_groups_df(adata_result, group="MG precursor").names[:10])
markers1 = list(sc.get.rank_genes_groups_df(adata_result, group="Neurogenesis").names[:10])
markers4 = list(sc.get.rank_genes_groups_df(adata_result, group="NRPC").names[:10])

time = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_10_monocle3_pseudotime.csv",sep = " ")
time["bin"] = pd.cut(time.x,2,labels = ["Early","Late"])

adata_result.obs['temp'] = adata_result.obs['temp'].astype(str)

adata_result.obs.loc[
    time.loc[time.bin=="Early"].index, "temp"
] = "Early Neurogenesis fate"
adata_result.obs.loc[
    time.loc[time.bin=="Late"].index, "temp"
] = "Late Neurogenesis fate"

#df = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/data/Human_TF_Full_List/TFsDatabaseExtract_v_1.01.csv")
#TF = df.loc[df["Is TF?"]=="Yes","HGNC symbol"]
#TF = [string for string in TF if not string.startswith("ZNF")]
#adata_result = adata_result[:,adata_result.var.index.isin(TF)]
temp = adata_result.copy()
#temp.obs.loc[temp.obs.temp=="MG precursor","temp"] = "Rest of RPC"
temp.obs["temp"] = temp.obs["temp"].astype(str)
temp2 = adata_result.copy()
temp2 = adata_result[adata_result.obs.temp.isin(["Early Neurogenesis fate","Late Neurogenesis fate"])]
sc.tl.rank_genes_groups(temp2, "temp")
markers2 = list(sc.get.rank_genes_groups_df(temp2, group="Early Neurogenesis fate").names[:10])
markers3 = list(sc.get.rank_genes_groups_df(temp2, group="Late Neurogenesis fate").names[:10])
plt.clf()
sc.pl.dotplot(
    temp,
    var_names={
        #"MG Precursor Enriched Genes": markers0,
        "Neurogenesis Fate Enriched Genes": markers1,
        "Early Neurogenesis Fate Enriched Genes": markers2,
        "Late Neurogenesis Fate Enriched Genes": markers3,
        "NRPC": markers4,
    },
    groupby="temp",
    dendrogram=False,
)
plt.xticks(rotation=90)
plt.title("")
plt.xlabel("")
plt.ylabel("")
fig = plt.gcf()
fig.set_size_inches(20, 5)
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

NRPC = adata_result[adata_result.obs.temp.isin(["Early Neurogenesis fate","Late Neurogenesis fate"]),:]
adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_MG_NRPC.h5ad")
NRPC = adata[NRPC.obs.index]
NRPC.obs["time"]=time.loc[NRPC.obs.index,"x"]
sc.pl.umap(NRPC,color = "time",frameon=False,size =50,vmax =2)
plt.title("")
plt.xlabel("")
plt.ylabel("")
fig = plt.gcf()
fig.set_size_inches(7, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/NRPC_fate_pseudo_time.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
    backend="cairo",
)
plt.clf()