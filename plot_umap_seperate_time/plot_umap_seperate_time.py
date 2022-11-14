import scvelo as scv
import plotly.express as px
import sys

input_adata = sys.argv[1]
output_fig_path = sys.argv[2]


def plot_adata(adata,title):
    df = adata.obs.copy()
    df["cell_label"] = adata.obs.index.values
    df["x"] = adata.obsm["X_umap"][:, 0]
    df["y"] = adata.obsm["X_umap"][:, 1]
    fig = px.scatter(
        df,
        x="x",
        y="y",
        title = title,
        color="scpred_prediction_mode",
        width=1000,
        height=1000,
        color_discrete_map={
            "AC": "#B6E880",
            "BC": "#EF553B",
            "Cone": "#00CC96",
            "HC": "#AB63FA",
            "MG": "#FFA15A",
            "RPC": "#19D3F3",
            "RGC": "#FF6692",
            "Rod": "#636EFA",
        },
    )
    fig.update_traces(mode="markers", marker_size=5)
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
    fig.update_layout(showlegend=False)
    fig.update_layout(
        margin=dict(l=0, r=0, b=0),
    )
    fig.update_layout(
    font=dict(
        size=30,
    )
)
    return fig


titles = ["FW10", "FW14", "FW16", "FW23"]
windows = [["10w", "11w2d"], ["14w2d", "14w5d"], ["16w4d","17w1d", "19w4d"], ["23w1d", "23w4d"]]
for i, window in enumerate(windows):
    adata = scv.read(input_adata)
    adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(str)
    if adata.obs.majorclass.str.contains("MG").any():
        adata.obs.at[adata.obs.majorclass == "MG", "scpred_prediction_mode"] = "MG"
    adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(
    "category"
)
    adata = adata[adata.obs.Time.isin(window)]
    fig = plot_adata(adata,titles[i])
    fig.write_html(output_fig_path + "cell_type_umap_" + str(i) + ".html")
    fig.write_image(output_fig_path + "cell_type_umap_" + str(i) + ".svg")
for i, window in enumerate(windows):
    adata = scv.read(input_adata)
    adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(str)
    if adata.obs.majorclass.str.contains("MG").any():
        adata.obs.at[adata.obs.majorclass == "MG", "scpred_prediction_mode"] = "MG"
    adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(
    "category"
)
    adata = adata[adata.obs.Region == "Macula"]
    adata = adata[adata.obs.Time.isin(window)]
    fig = plot_adata(adata,"Macula " + titles[i])
    fig.write_html(output_fig_path + "cell_type_umap_" + str(i) + "_Macula.html")
    fig.write_image(output_fig_path + "cell_type_umap_" + str(i) + "_Macula.svg")
    
for i, window in enumerate(windows):
    adata = scv.read(input_adata)
    adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(str)
    if adata.obs.majorclass.str.contains("MG").any():
        adata.obs.at[adata.obs.majorclass == "MG", "scpred_prediction_mode"] = "MG"
    adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(
    "category"
)
    adata = adata[adata.obs.Region == "Peripheral"]
    adata = adata[adata.obs.Time.isin(window)]
    fig = plot_adata(adata,"Peripheral " + titles[i])
    fig.write_html(output_fig_path + "cell_type_umap_" + str(i) + "_Peripheral.html")
    fig.write_image(output_fig_path + "cell_type_umap_" + str(i) + "_Peripheral.svg")
