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
from scipy.cluster.hierarchy import leaves_list


hs = joblib.load(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/PPRC_hs.create_modules.pkl"
)
adata = joblib.load(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/PPRC_hs.adata.pkl"
)


modules = hs.create_modules(min_gene_threshold=200, core_only=True, fdr_threshold=0.05)
modules.to_csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/modules.csv")

def my_plot_local_correlations(
        self, mod_cmap="tab10", vmin=-8, vmax=8, z_cmap="RdBu_r", yticklabels=False
    ):
        """Plots a clustergrid of the local correlation values
        Parameters
        ----------
        mod_cmap: valid matplotlib colormap str or object
            discrete colormap for module assignments on the left side
        vmin: float
            minimum value for colorscale for Z-scores
        vmax: float
            maximum value for colorscale for Z-scores
        z_cmap: valid matplotlib colormap str or object
            continuous colormap for correlation Z-scores
        yticklabels: bool
            Whether or not to plot all gene labels
            Default is false as there are too many.  However
            if using this plot interactively you may with to set
            to true so you can zoom in and read gene names
        """

        return my_local_correlation_plot(
            self.local_correlation_z,
            self.modules,
            self.linkage,
            mod_cmap=mod_cmap,
            vmin=vmin,
            vmax=vmax,
            z_cmap=z_cmap,
            yticklabels=yticklabels,
        )

def my_local_correlation_plot(
            local_correlation_z, modules, linkage,
            mod_cmap='tab10', vmin=-8, vmax=8,
            z_cmap='RdBu_r', yticklabels=True
):

    row_colors = None
    colors = list(plt.get_cmap(mod_cmap).colors)
    module_colors = {i: colors[(i-1) % len(colors)] for i in modules.unique()}
    module_colors[-1] = '#ffffff'

    row_colors1 = pd.Series(
        [module_colors[i] for i in modules],
        index=local_correlation_z.index,
    )

    row_colors = pd.DataFrame({
        "Modules": row_colors1,
    })

    cm = sns.clustermap(
        local_correlation_z,
        row_linkage=linkage,
        col_linkage=linkage,
        vmin=vmin,
        vmax=vmax,
        cmap=z_cmap,
        xticklabels=False,
        yticklabels=yticklabels,
        row_colors=row_colors,
        rasterized=True,
    )
    hm = cm.ax_heatmap
    y_labels = [tick.get_text() for tick in hm.get_yticklabels()]
    modules[y_labels].to_csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/hm.get_yticklabels.csv")
    fig = plt.gcf()
    plt.sca(cm.ax_heatmap)
    plt.ylabel("")
    plt.xlabel("")

    cm.ax_row_dendrogram.remove()

    # Add 'module X' annotations
    ii = leaves_list(linkage)

    mod_reordered = modules.iloc[ii]

    mod_map = {}
    y = np.arange(modules.size)

    for x in mod_reordered.unique():
        if x == -1:
            continue

        mod_map[x] = y[mod_reordered == x].mean()

    plt.sca(cm.ax_row_colors)
    for mod, mod_y in mod_map.items():
        plt.text(-.5, y=mod_y, s="Module {}".format(mod),
                 horizontalalignment='right',
                 verticalalignment='center')
    plt.xticks([])

    # Find the colorbar 'child' and modify
    min_delta = 1e99
    min_aa = None
    for aa in fig.get_children():
        try:
            bbox = aa.get_position()
            delta = (0-bbox.xmin)**2 + (1-bbox.ymax)**2
            if delta < min_delta:
                delta = min_delta
                min_aa = aa
        except AttributeError:
            pass

    min_aa.set_ylabel('Z-Scores')
    min_aa.yaxis.set_label_position("left")

my_local_correlation_plot(hs.local_correlation_z,
            hs.modules,
            hs.linkage)

hs.plot_local_correlations()
plt.savefig(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/plot_local_correlations.svg"
)


module_scores = hs.calculate_module_scores()

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_run_umap_PRPC/adata_umap.h5ad"
)

adata_result.obs["Module1"] = module_scores.loc[adata_result.obs.index,1].values
adata_result.obs["Module2"] = module_scores.loc[adata_result.obs.index,2].values
adata_result.obs["Module3"] = module_scores.loc[adata_result.obs.index,3].values
adata_result.obs["Module4"] = module_scores.loc[adata_result.obs.index,4].values

adata_result.write(
    "/storage/singlecell/zz4/fetal_bash/temp/adata_umap.h5ad"
)

module_list = [
    "Module1",
    "Module2",
    "Module3",
    "Module4",
]
cols = [
    "#636EFA",
    "#EF553B",
    "#00CC96",
    "#FFA15A",
    "#AB63FA",
    "#19D3F3",
    "#FF6692",
    "#B6E880",
    "#FF97FF",
    "#FECB52",
]
for i in range(len(module_list)):
    width = 1000
    height = 1000
    legend_size = 3
    marker_size = 8
    df = adata_result.obs.copy()
    df["cell_label"] = adata_result.obs.index.values
    df["x"] = adata_result.obsm["X_umap"][:, 0]
    df["y"] = adata_result.obsm["X_umap"][:, 1]
    fig = px.scatter(
        df,
        x="x",
        y="y",
        color=module_list[i],
        width=width,
        height=height,
        color_continuous_scale=["#BAB0AC", cols[i]],
        hover_data=["cell_label"],
        range_color=[0.5, 1.5],
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
    fig.write_image(
        "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/PRPC_gene_"
        + module_list[i]
        + ".svg"
    )
