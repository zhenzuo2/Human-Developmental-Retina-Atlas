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
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_cell_type.svg",
    dpi=600,
)

scv.pl.umap(
    adata,
    color="Region",
    size=1,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_region.svg",
    dpi=600,
)

adata.obs["Days"] = adata.obs["Days"].astype(float)
scv.pl.umap(
    adata,
    color="Days",
    size=1,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_days.svg",
    dpi=600,
)
