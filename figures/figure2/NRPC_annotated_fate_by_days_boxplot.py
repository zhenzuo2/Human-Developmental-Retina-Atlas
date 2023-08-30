# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns

sc.set_figure_params(transparent=True, fontsize=40)
NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC = NRPC[NRPC.obs.leiden != "12"]
NRPC.obs["subclass"] = 'Undetermined'

sc.pp.neighbors(NRPC, use_rep="X_scVI")
sc.tl.leiden(NRPC, resolution=3)

clusters = ["28", "3", "13", "15"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["6", "30", "35", "25"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["21", "14", "27"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["33", "0"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["12", "22"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["4", "19", "23"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"

NRPC.obs["Days"] = NRPC.obs["Days"].astype(float)
NRPC.obs["PCW"] = NRPC.obs["Days"] / 7

# Create a boxplot using Seaborn
sns.boxplot(
    data=NRPC.obs,
    x="subclass",
    y="PCW",
    hue="Region",
    order=['Undetermined',"RGC", "Cone", "HC", "AC", "Rod", "BC"],
    showfliers=False,
)
plt.title("")
plt.xlabel("Precursor Group")
plt.ylabel("PCW")
plt.xticks(rotation = 90)
# Move the legend outside the figure
plt.legend(title="Region", loc="best", bbox_to_anchor=(1, 0.5))
fig = plt.gcf()
fig.set_size_inches(15, 15)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_subclass_Days.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()