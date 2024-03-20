import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

NRPC = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC = NRPC[NRPC.obs.leiden != "12"]
NRPC.obs["subclass"] = 'Undetermined'

sc.pp.neighbors(NRPC, use_rep="X_scVI")
sc.tl.leiden(NRPC, resolution=3)

clusters = ["28", "3", "13", "15"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"

adata_result = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad")
adata_result = adata_result[NRPC[NRPC.obs.subclass=="AC"].obs.index,]
sc.pp.normalize_total(adata_result)
sc.pp.log1p(adata_result)


adata_result.obs["Weeks"] = adata_result.obs.Days.map(
    {
        70: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW16",
        103: "PCW16",
        116: "PCW16",
        120: "PCW16",
        136: "PCW20",
        137: "PCW20",
        141: "PCW20",
        142: "PCW20",
        162: "PCW23",
        165: "PCW23",
    }
)

sc.pl.violin(adata_result[adata_result.obs.Weeks.isin(['PCW10', 'PCW13', 'PCW16', 'PCW20'])], keys='ONECUT2',groupby = "Weeks")
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_ONECUT2.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    #backend="cairo",
)
