# Import packages
import scvelo as scv
import matplotlib.pyplot as plt

# set verbosity levels
cr.settings.verbosity = 2
scv.settings.verbosity = 3

scv.set_figure_params(dpi=200, dpi_save=600)

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000.h5ad"
)
adata.obsm["X_umap"]
umap = pd.DataFrame(adata.obsm["X_umap"], columns=["x", "y"])
umap.index = adata.obs.index
umap

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)

adata.obsm["X_umap"] = umap.loc[
    adata.obs.index,
].values

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
plt.savefig("/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_cell_type.svg")