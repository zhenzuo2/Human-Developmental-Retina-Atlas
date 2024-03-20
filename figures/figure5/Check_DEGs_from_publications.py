import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib

matplotlib.rcParams.update({"font.size": 18})
genes = [
    "TYR",
    "OCA2",
    "TYRP1",
    "SLC45A2",
    "SLC24A5",
    "LRMDA",
    "GPR143",
    "HPS1",
    "AP3B1",
    "HPS3",
    "HPS4",
    "HPS5",
    "HPS6",
    "DTNBP1",
    "BLOC1S3",
    "BLOC1S6",
    "AP3D1",
    "BLOC1S5",
    "LYST",
    "SLC38A8",
    "PAX6",
    "FRMD7",
    "AHR",
    "CNGB3",
    "CNGA3",
    "GNAT2",
    "PDE6C",
    "PDE6H",
    "ATF6",
]

adata = sc.read("merged_raw_filtered_umap_10000_woadult_MG.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
for x in ["AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod"]:
    print(x)
    temp = adata[adata.obs.majorclass == x]
    sc.pp.highly_variable_genes(temp, n_top_genes=5000, subset=True)
    sc.tl.rank_genes_groups(temp, groupby="Region")
    df = sc.get.rank_genes_groups_df(temp, group="Macula")
    df.index = df.names.values
    print(
        df.loc[[x for x in genes if x in temp.var.index],]
        .loc[df.pvals_adj < 0.01,]
        .sort_values("scores")
    )

x = "Cone"
print(x)
temp = adata[adata.obs.majorclass == x]
sc.pp.highly_variable_genes(temp, n_top_genes=2000, subset=True)
sc.tl.rank_genes_groups(temp, groupby="Region")
df = sc.get.rank_genes_groups_df(temp, group="Macula")
df.index = df.names.values
print(
    df.loc[[x for x in genes if x in temp.var.index],]
    .loc[df.pvals_adj < 0.01,]
    .sort_values("scores")
)
M_cells = (
    adata[(adata.obs.majorclass == x) & (adata.obs.Region == "Macula")]
    .obs.sample(n=2000)
    .index
)
P_cells = (
    adata[(adata.obs.majorclass == x) & (adata.obs.Region == "Peripheral")]
    .obs.sample(n=2000)
    .index
)
adata = adata[list(M_cells) + list(P_cells),]
df = sc.get.obs_df(
    adata[adata.obs.majorclass == x],
    list(
        df.loc[[x for x in genes if x in temp.var.index],]
        .loc[df.pvals_adj < 0.01,]
        .sort_values("scores")
        .names
    )
    + ["Region"],
)
df["Region"] = adata[adata.obs.majorclass == x].obs.Region
df = df.set_index("Region").stack().reset_index()
df.columns = ["Region", "gene", "value"]
import seaborn as sns

sns.violinplot(data=df, x="gene", y="value", hue="Region", split=True)
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
sns.swarmplot(
    x="gene",
    y="value",
    data=df.loc[df.value > 0],
    hue="Region",
    size=2,
    dodge=True,
    legend=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "Cone.png",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()
