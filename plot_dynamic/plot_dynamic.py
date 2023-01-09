import scvelo as scv
import os

try:
    os.makedirs(
        os.path.dirname(
            "/storage/singlecell/zz4/fetal_bash/figures/scv.pl.velocity_embedding_stream/"
        )
    )
except FileExistsError:
    # directory already exists
    pass

cell_types = [
    "AC",
    "all",
    "BC",
    "Cone",
    "HC",
    "MG",
    "NRPC",
    "PRPC",
    "RGC",
    "Rod",
    "RPC",
]
color_list = ["scpred_prediction", "Time", "Region", "Days", "subclass", "majorclass"]
for cell_type in cell_types:
    adata = scv.read(
        "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/"
        + cell_type
        + ".h5ad"
    )
    adata.obs['Days'] = adata.obs['Days'].astype(float)
    for color in color_list:
        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            legend_fontsize=20,
            title="",
            color=color,
            size=20,
            colorbar=True,
            dpi=600,
            figsize=(10, 10),
            linewidth=2,
            alpha=1,
            legend_loc="on data",
            save="/storage/singlecell/zz4/fetal_bash/figures/scv.pl.velocity_embedding_stream/"
            + cell_type
            + "_"
            + color
            + "_scv.pl.velocity_embedding_stream.png",
        )

