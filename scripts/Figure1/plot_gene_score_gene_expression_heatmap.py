# Import packages
import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import magic


def smooth(adata):
    magic_op = magic.MAGIC()
    magic_op.set_params(n_jobs=10)
    emt_magic = magic_op.fit_transform(adata.X, genes="all_genes")
    emt_magic = magic_op.transform(genes="all_genes")
    adata.X = emt_magic
    return adata


plt.rcParams["font.family"] = "Arial"
scv.set_figure_params(dpi=300, dpi_save=300)
sc.set_figure_params(scanpy=True, dpi_save=300, fontsize=25)

gs = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score.h5ad"
)
adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=True)

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

sc.pp.normalize_total(gs, target_sum=1e6)
sc.pp.log1p(gs)

gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]

cells = [x for x in adata.obs.index if x in gs.obs.index]
adata = adata[cells]
gs = gs[cells]
adata.obs = pd.DataFrame(adata.obs)

gs.obs = adata.obs

sc.tl.rank_genes_groups(adata, "majorclass")

output = "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table3/Supplementary Table 3.xlsx"
df = sc.get.rank_genes_groups_df(adata, group="BC")
with pd.ExcelWriter(output) as writer:
    df.to_excel(writer, sheet_name="BC", index=False)
BC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
BC = [x for x in BC if x in gs.var.index][:50]
random.shuffle(BC)

df = sc.get.rank_genes_groups_df(adata, group="Cone")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="Cone", index=False)
Cone = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
Cone = [x for x in Cone if x in gs.var.index][:50]
random.shuffle(Cone)

df = sc.get.rank_genes_groups_df(adata, group="Rod")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="Rod", index=False)
Rod = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
Rod = [x for x in Rod if x in gs.var.index][:50]
random.shuffle(Rod)

df = sc.get.rank_genes_groups_df(adata, group="MG")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="MG", index=False)
MG = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
MG = [x for x in MG if x in gs.var.index][:50]
random.shuffle(MG)

df = sc.get.rank_genes_groups_df(adata, group="PRPC")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="PRPC", index=False)
PRPC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
PRPC = [x for x in PRPC if x in gs.var.index][:50]
random.shuffle(PRPC)

df = sc.get.rank_genes_groups_df(adata, group="NRPC")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="NRPC", index=False)
NRPC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
NRPC = [x for x in NRPC if x in gs.var.index][:50]
random.shuffle(NRPC)

df = sc.get.rank_genes_groups_df(adata, group="RGC")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="RGC", index=False)
RGC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
RGC = [x for x in RGC if x in gs.var.index][:50]
random.shuffle(RGC)

df = sc.get.rank_genes_groups_df(adata, group="AC")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="AC", index=False)
AC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
AC = [x for x in AC if x in gs.var.index][:50]
random.shuffle(AC)

df = sc.get.rank_genes_groups_df(adata, group="HC")
with pd.ExcelWriter(output, engine="openpyxl", mode="a") as writer:
    df.to_excel(writer, sheet_name="HC", index=False)
HC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
HC = [x for x in HC if x in gs.var.index][:50]
random.shuffle(HC)

Markers = {
    "BC": BC,
    "Cone": Cone,
    "Rod": Rod,
    "MG": MG,
    "NRPC": NRPC,
    "PRPC": PRPC,
    "RGC": RGC,
    "AC": AC,
    "HC": HC,
}
df = pd.DataFrame.from_dict(Markers, orient="index").transpose()
df.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/Top50_Markers.csv")

nested_list = list(Markers.values())
flattened_list = list(set([item for sublist in nested_list for item in sublist]))

df = pd.DataFrame(gs.obs)
df.loc[:, "cell_id"] = list(df.index.values)
grouped_data = df.groupby("majorclass")
# Define the number of rows to downsample
num_rows = 4000
# Sample the same number of rows from each group
downsampled_data = grouped_data.apply(
    lambda x: x.sample(n=num_rows, random_state=42, replace=True)
)
# Reset the index
downsampled_data = downsampled_data.reset_index(drop=True)

gs_subset = gs[downsampled_data.cell_id]
gs_subset.obs["majorclass"] = pd.Categorical(
    list(gs_subset.obs["majorclass"]),
    categories=["BC", "Cone", "Rod", "MG", "NRPC", "PRPC", "RGC", "AC", "HC"],
)
gs_subset = gs_subset[:, flattened_list]
sc.pp.scale(gs_subset)
gs_subset = smooth(gs_subset)
gs_subset.write(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/gene_score_heatmap.h5ad"
)
sc.pl.heatmap(
    gs_subset,
    Markers,
    groupby="majorclass",
    cmap="plasma",
    dendrogram=False,
    show_gene_labels=False,
    vmax=0.65,
    vmin=0.55,
    standard_scale="obs",
)
fig = plt.gcf()
fig.set_size_inches(12, 10)
plt.xlabel("")
plt.ylabel("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/gene_score_heatmap.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
##########################################################################################################################
adata_subset = adata[downsampled_data.cell_id]
adata_subset = adata_subset[:, flattened_list]
sc.pp.scale(adata_subset)
adata_subset = smooth(adata_subset)
adata_subset.obs["majorclass"] = pd.Categorical(
    list(adata_subset.obs["majorclass"]),
    categories=["BC", "Cone", "Rod", "MG", "NRPC", "PRPC", "RGC", "AC", "HC"],
)
sc.pl.heatmap(
    adata_subset,
    Markers,
    groupby="majorclass",
    dendrogram=False,
    show_gene_labels=False,
    vmax=0.6,
    vmin=0.1,
    cmap="viridis",
    standard_scale="obs",
)
fig = plt.gcf()
plt.xlabel("")
plt.ylabel("")
fig.set_size_inches(12, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/gene_expression_heatmap.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)
adata_subset.write(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/gene_expression_heatmap.h5ad"
)
