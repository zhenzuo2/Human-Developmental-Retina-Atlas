import sys
import scvelo as scv
import plotly.express as px
import os

adata_file = sys.argv[1]
output_fig_path = sys.argv[2]

# Check whether the specified path exists or not
isExist = os.path.exists(output_fig_path)

if not isExist:
  
  # Create a new directory because it does not exist 
  os.makedirs(output_fig_path)

adata = scv.read(adata_file)

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
    fig.write_image(output_fig_path + attribute + ".svg")