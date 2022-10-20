import sys
import scvelo as scv
import pandas as pd
import plotly.express as px

adata_file = sys.argv[1]
cell_type = sys.argv[2]
meta_cluster_adata_file = sys.argv[3]
subclass_reference_file = sys.argv[4]
majorclass_reference_file = sys.argv[5]
output_file_path = sys.argv[6]
output_fig_path = sys.argv[7]

adata = scv.read(adata_file)

adata = adata[
    adata.obs.scpred_prediction_mode == cell_type,
]

# Remove this unwanted cluster of cells
if cell_type == "AC":
    adata = adata[adata.obs.leiden != "16"]

if cell_type == "BC":
    adata = adata[adata.obs.leiden != "28"]
    adata = adata[adata.obs.leiden != "41"]

if cell_type == "RGC":
    adata = adata[adata.obs.leiden != "8"]
    adata = adata[adata.obs.leiden != "29"]

# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs.loc[adata.obs.index, :]
meta_cluster["author_cell_type"] = meta_cluster.author_cell_type.astype(str).astype(int)

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.map(
    dict(zip(subclass_reference.leiden, subclass_reference.subclass))
).fillna(cell_type + " Precursor")

# Read reference file
majorclass_reference = pd.read_csv(majorclass_reference_file)
adata.obs["majorclass"] = (
    adata.obs["subclass"]
    .map(dict(zip(majorclass_reference.subclass, majorclass_reference.majorclass)))
    .fillna(cell_type + " Precursor")
)

adata.write(output_file_path + cell_type + "_major_sub_class.h5ad")
adata.obs.to_csv(output_file_path + cell_type + "_major_sub_class_obs.csv")

for attribute in ["Days", "Region", "subclass", "majorclass"]:
    df = adata.obs.copy()
    df["x"] = adata.obsm["X_umap"][:, 0]
    df["y"] = adata.obsm["X_umap"][:, 1]
    fig = px.scatter(
        df,
        x="x",
        y="y",
        color=attribute,
        width=1200,
        height=1000,
    )
    fig.update_traces(mode="markers", marker_size=3)
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)
    fig.update_layout(yaxis_visible=False, yaxis_showticklabels=False)
    fig.update_layout(xaxis_visible=False, xaxis_showticklabels=False)
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0, 0, 0, 0)",
            "paper_bgcolor": "rgba(0, 0, 0, 0)",
        }
    )
    fig.update_layout(legend={"title_text": ""})
    fig.update_layout(legend={"itemsizing": "constant"})
    fig.update_layout(legend=dict(font=dict(size=5)))
    fig.write_image(output_fig_path + cell_type + attribute + ".svg")