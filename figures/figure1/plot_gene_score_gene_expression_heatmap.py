# Import packages
import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

scv.set_figure_params(dpi=600, dpi_save=600)

gs = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/gene_score.h5ad"
)
adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad"
)
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
sc.pp.scale(adata)
sc.pp.scale(gs)

sc.tl.rank_genes_groups(adata, 'majorclass')

BC = (
    sc.get.rank_genes_groups_df(adata, group="BC")
    .loc[sc.get.rank_genes_groups_df(adata, group="BC").pvals_adj < 0.05, "names"]
    .values
)
BC = [x for x in BC if x in gs.var.index]

Cone = (
    sc.get.rank_genes_groups_df(adata, group="Cone")
    .loc[sc.get.rank_genes_groups_df(adata, group="Cone").pvals_adj < 0.05, "names"]
    .values
)
Cone = [x for x in Cone if x in gs.var.index]

Rod = (
    sc.get.rank_genes_groups_df(adata, group="Rod")
    .loc[sc.get.rank_genes_groups_df(adata, group="Rod").pvals_adj < 0.05, "names"]
    .values
)
Rod = [x for x in Rod if x in gs.var.index]

MG = (
    sc.get.rank_genes_groups_df(adata, group="MG")
    .loc[sc.get.rank_genes_groups_df(adata, group="MG").pvals_adj < 0.05, "names"]
    .values
)
MG = [x for x in MG if x in gs.var.index]

RPC = (
    sc.get.rank_genes_groups_df(adata, group="RPC")
    .loc[sc.get.rank_genes_groups_df(adata, group="RPC").pvals_adj < 0.05, "names"]
    .values
)
RPC = [x for x in RPC if x in gs.var.index]

RGC = (
    sc.get.rank_genes_groups_df(adata, group="RGC")
    .loc[sc.get.rank_genes_groups_df(adata, group="RGC").pvals_adj < 0.05, "names"]
    .values
)
RGC = [x for x in RGC if x in gs.var.index]

AC = (
    sc.get.rank_genes_groups_df(adata, group="AC")
    .loc[sc.get.rank_genes_groups_df(adata, group="AC").pvals_adj < 0.05, "names"]
    .values
)
AC = [x for x in AC if x in gs.var.index]

HC = (
    sc.get.rank_genes_groups_df(adata, group="HC")
    .loc[sc.get.rank_genes_groups_df(adata, group="HC").pvals_adj < 0.05, "names"]
    .values
)
HC = [x for x in HC if x in gs.var.index]

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
num_rows = 100
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
    vmax=np.quantile(gs_subset.X, 0.95),
)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/gene_score_heatmap.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)
##########################################################################################################################
adata_subset = adata[downsampled_data.cell_id]
ax = sc.pl.heatmap(
    adata_subset,
    Markers,
    groupby="majorclass",
    dendrogram=False,
    show_gene_labels=False,
    vmax=np.quantile(adata_subset.X.toarray(), 0.99),
    cmap="plasma",
)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/gene_expression_heatmap.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
    backend="cairo",
)
