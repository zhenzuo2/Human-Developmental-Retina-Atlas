import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/AC.h5ad")
adata_result = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad")

meta_file = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/AC_subtype.csv")
meta_file.index = meta_file["Unnamed: 0"]
common_cells = [x for x in meta_file.index if x in adata.obs.index]
meta_file = meta_file.loc[common_cells, :]
adata.obs["subclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
for sublcass in set(meta_file.subclass):
    adata.obs.loc[
        meta_file.loc[meta_file.subclass == sublcass, :].index, "subclass"
    ] = sublcass
for majorclass in set(meta_file.majorclass):
    adata.obs.loc[
        meta_file.loc[meta_file.majorclass == majorclass, :].index, "majorclass"
    ] = majorclass

adata_result = adata_result[adata.obs.index]
adata_result.obsm["X_umap"] = adata.obsm["X_umap"]
adata_result.obs["leiden"] = adata.obs["leiden"]

sc.pp.normalize_per_cell(adata_result)
sc.pp.log1p(adata_result)

precursors = adata_result[adata_result.obs.leiden.isin(["7", "5", "8"])]
precursors.obs["leiden"] = precursors.obs.leiden.replace(
    {"7": "Starburst", "5": "GABAergic", "8": "Glycinergic"}
)
sc.tl.rank_genes_groups(precursors, "leiden")

sc.tl.dendrogram(precursors, groupby="leiden")
sc.pl.rank_genes_groups_heatmap(
    precursors,
    n_genes=20,
    use_raw=False,
    swap_axes=True,
    vmin=-3,
    vmax=3,
    cmap="bwr",
    figsize=(10, 10),
    show=False,
    show_gene_labels=True,
)
plt.xticks(color="w")
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_precursor_de_gene.svg",
    dpi=600,
    bbox_inches="tight",
)
