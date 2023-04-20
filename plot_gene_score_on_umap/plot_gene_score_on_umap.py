import scvelo as scv
import scanpy as sc
gs = scv.read("/storage/singlecell/zz4/fetal_bash/results/Gene_score_matrix/gene_score.h5ad")
adata_result = scv.read("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad")

gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]

cells = [x for x in adata_result.obs.index if x in gs.obs.index]
adata_result = adata_result[cells]
gs = gs[cells]

gs.obsm["X_umap"] = adata_result.obsm["X_umap"]

scv.pp.normalize_per_cell(adata_result)
scv.pp.log1p(adata_result)

sc.pl.umap(
    adata_result,
    color="NR2E3",
    #vmax=3,
    frameon=False,
    size=3,
    save="NR2E3.svg",
    title="",
)

sc.pl.umap(
    gs,
    color="NR2E3",
    vmax=3,
    frameon=False,
    size=3,
    save="NR2E3.svg",
    title="",
)