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

color_list = ["scpred_prediction", "Time", "Region", "Days", "subclass", "majorclass"]
adata = scv.read("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class_ldata_dynamic.h5ad")
adata.obs['Days'] = adata.obs['Days'].astype(float)
adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "AC Precursor": "AC",
        "BC Precursor": "BC",
        "Cone Precursor": "Cone",
        "GABAergic": "AC",
        "Glycinergic": "AC",
        "HC0": "HC",
        "HC1": "HC",
        "MG": "MG",
        "ML_Cone": "Cone",
        "NRPC": "RPC",
        "OFF-BC": "BC",
        "OFF_MGC": "RGC",
        "ON-BC": "BC",
        "ON_MGC": "RGC",
        "PRPC": "RPC",
        "RBC": "BC",
        "RGC Precursor": "RGC",
        "Rod": "Rod",
        "Rod Precursor": "Rod",
        "SACs": "AC",
        "S_Cone": "Cone",
        "dual ACs": "AC",
    }
)

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
        + color
        + "_scv.pl.velocity_embedding_stream.png",
    )

