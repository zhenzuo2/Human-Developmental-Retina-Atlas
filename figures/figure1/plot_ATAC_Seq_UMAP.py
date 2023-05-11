# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)

adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)
adata.obs["scpred_prediction"] = pd.Categorical(
    list(adata.obs["scpred_prediction"]),
    categories=[
        "RPC",
        "RGC",
        "Cone",
        "HC",
        "AC",
        "Rod",
        "BC",
        "MG",
    ],
)

meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/atac_umap.csv", sep=" "
)
common_cells = [x for x in meta.index if x in adata.obs.index]
adata = adata[common_cells]
meta = meta.loc[common_cells, :]
meta["IterativeLSI#UMAP_Dimension_2"] = -meta["IterativeLSI#UMAP_Dimension_2"]
adata.obsm["X_umap"] = meta[
    ["IterativeLSI#UMAP_Dimension_1", "IterativeLSI#UMAP_Dimension_2"]
].to_numpy()

scv.pl.umap(
    adata,
    color="scpred_prediction",
    size=1,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_atac_umap_by_cell_type.svg",
    dpi=600,
    bbox_inches="tight",
)
