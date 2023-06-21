import joblib
import scvelo as scv
import plotly.express as px
import scanpy as sc
import hotspot
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib

matplotlib.rcParams.update({"font.size": 30})

output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"
hs = joblib.load(
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/latent_time/PPRC_hs.create_modules.pkl"
)
adata = joblib.load(
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/latent_time/PPRC_hs.adata.pkl"
)

modules = hs.create_modules(min_gene_threshold=100, core_only=True, fdr_threshold=0.05)

hs.modules = hs.modules.replace({4: 1, 5: 2, 1: 4, 2: 5})
hs.modules.to_csv(output_file_path + "PRPC_gene_modules.csv")
hs.plot_local_correlations()
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "PRPC_gene_module_heatmap.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
)

module_scores = hs.calculate_module_scores()
adata_result = joblib.load(
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/latent_time/PPRC_hs.adata.pkl"
)
PRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC/adata_umap.h5ad"
)
adata_result.obsm["X_umap"] = PRPC[adata_result.obs.index].obsm["X_umap"]

adata_result.obs["Module1"] = module_scores.loc[adata_result.obs.index, 1].values
adata_result.obs["Module2"] = module_scores.loc[adata_result.obs.index, 2].values
adata_result.obs["Module3"] = module_scores.loc[adata_result.obs.index, 3].values
adata_result.obs["Module4"] = module_scores.loc[adata_result.obs.index, 4].values
adata_result.obs["Module5"] = module_scores.loc[adata_result.obs.index, 5].values
adata_result.obs["Module"] = adata_result.obs[
    ["Module1", "Module2", "Module3", "Module4", "Module5"]
].idxmax(axis=1)

module_list = [
    "Module1",
    "Module2",
    "Module3",
    "Module4",
    "Module5",
]
cols = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
]
for i in range(len(module_list)):
    width = 1200
    height = 1200
    legend_size = 3
    marker_size = 7
    adata_result.obs["temp"] = adata_result.obs.Module == module_list[i]
    adata_result.obs["temp"] = adata_result.obs["temp"].replace(
        {True: "True", False: "False"}
    )
    adata_result.obs["temp"].astype("category").cat.reorder_categories(
        ["True", "False"], inplace=True
    )
    df = adata_result.obs.copy()
    df["cell_label"] = adata_result.obs.index.values
    df["x"] = list(adata_result.obsm["X_umap"][:, 0])
    df["y"] = list(adata_result.obsm["X_umap"][:, 1])
    df = df.sort_values(by=[module_list[i]])
    fig = px.scatter(
        df,
        x="x",
        y="y",
        color="temp",
        width=width,
        height=height,
        color_discrete_map={"True": cols[i], "False": "lightgray"},
        hover_data=["cell_label"],
        range_color=[0.5, 4],
    )
    fig.update_traces(mode="markers", marker_size=marker_size)
    fig.update_layout(showlegend=False)
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
    fig.update_layout(legend=dict(font=dict(size=legend_size)))
    fig.write_image(output_file_path + module_list[i] + ".svg")
