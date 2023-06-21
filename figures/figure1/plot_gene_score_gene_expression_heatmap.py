# Import packages
import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
scv.set_figure_params(dpi=600, dpi_save=600)

gs = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/gene_score.h5ad"
)
adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad"
)
sc.pp.filter_genes(adata, min_counts=10000)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=True)

scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)

gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]

cells = [x for x in adata.obs.index if x in gs.obs.index]
adata = adata[cells]
gs = gs[cells]
adata.obs = pd.DataFrame(adata.obs)
adata.obs["majorclass"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)
gs.obs = adata.obs

sc.tl.rank_genes_groups(adata, "majorclass",method = "wilcoxon")

sc.pp.scale(adata)
sc.pp.scale(gs)

df = sc.get.rank_genes_groups_df(adata, group="BC")
BC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
BC = [x for x in BC if x in gs.var.index]
random.shuffle(BC)

df = sc.get.rank_genes_groups_df(adata, group="Cone")
Cone = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
Cone = [x for x in Cone if x in gs.var.index]
random.shuffle(Cone)

df = sc.get.rank_genes_groups_df(adata, group="Rod")
Rod = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
Rod = [x for x in Rod if x in gs.var.index]
random.shuffle(Rod)

df = sc.get.rank_genes_groups_df(adata, group="MG")
MG = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
MG = [x for x in MG if x in gs.var.index]
random.shuffle(MG)

df = sc.get.rank_genes_groups_df(adata, group="RPC")
RPC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
RPC = [x for x in RPC if x in gs.var.index]
random.shuffle(RPC)

df = sc.get.rank_genes_groups_df(adata, group="RGC")
RGC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
RGC = [x for x in RGC if x in gs.var.index]
random.shuffle(RGC)

df = sc.get.rank_genes_groups_df(adata, group="AC")
AC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
AC = [x for x in AC if x in gs.var.index]
random.shuffle(AC)

df = sc.get.rank_genes_groups_df(adata, group="HC")
HC = df.loc[(df.pvals_adj < 0.01) & (df.logfoldchanges >= 2), "names"].values
HC = [x for x in HC if x in gs.var.index]
random.shuffle(HC)

Markers = {
    "BC": BC,
    "Cone": Cone,
    "Rod": Rod,
    "MG": MG,
    "RPC": RPC,
    "RGC": RGC,
    "AC": AC,
    "HC": HC,
}

df = pd.DataFrame(gs.obs)
df.loc[:, "cell_id"] = list(df.index.values)
grouped_data = df.groupby("majorclass")
# Define the number of rows to downsample
num_rows = 10000
# Sample the same number of rows from each group
downsampled_data = grouped_data.apply(
    lambda x: x.sample(n=num_rows, random_state=42, replace=True)
)
# Reset the index
downsampled_data = downsampled_data.reset_index(drop=True)

sc.set_figure_params(scanpy=True, dpi_save=600, fontsize=25)
gs_subset = gs[downsampled_data.cell_id]
gs_subset.obs["majorclass"] = pd.Categorical(
    list(gs_subset.obs["majorclass"]),
    categories=["BC", "Cone", "Rod", "MG", "RPC", "RGC", "AC", "HC"],
)
sc.pl.heatmap(
    gs_subset,
    Markers,
    groupby="majorclass",
    cmap="viridis",
    dendrogram=False,
    show_gene_labels=False,
    vmax=np.quantile(gs_subset.X, 0.9),
    vmin=np.quantile(gs_subset.X, 0.1),
)
fig = plt.gcf()
fig.set_size_inches(12, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/gene_score_heatmap.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)
##########################################################################################################################
adata_subset = adata[downsampled_data.cell_id]
adata_subset.obs["majorclass"] = pd.Categorical(
    list(adata_subset.obs["majorclass"]),
    categories=["BC", "Cone", "Rod", "MG", "RPC", "RGC", "AC", "HC"],
)
ax = sc.pl.heatmap(
    adata_subset,
    Markers,
    groupby="majorclass",
    dendrogram=False,
    show_gene_labels=False,
    vmax=np.quantile(adata_subset.X, 0.95),
    vmin=np.quantile(adata_subset.X, 0.05),
    cmap="plasma",
)
fig = plt.gcf()
fig.set_size_inches(12, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/gene_expression_heatmap.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
    backend="cairo",
)
