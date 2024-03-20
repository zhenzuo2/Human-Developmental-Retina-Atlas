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
from scipy.cluster.hierarchy import leaves_list

matplotlib.rcParams.update({"font.size": 28})

output_file_path = "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/"
hs = joblib.load(
    "/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/PRPC_hs.create_modules_velocity_pseudotime.pkl"
)

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

modules = hs.create_modules(min_gene_threshold=160, core_only=True, fdr_threshold=0.01)
hs.modules = hs.modules.replace({3: 1, 4: 2, 2: 3})
for x in [
    "GNB4",
    "TEAD4",
    "HAUS2",
    "VKORC1L1",
    "LMNB1-DT",
    "PCGF6",
    "EIF2AK3",
    "FAM185A",
    "AC122719.3",
    "NEDD4",
    "LRRC8B",
    "SPDYA",
    "AC092910.3",
    "SV2C",
    "SOX11",
    "ABHD17B",
]:
    hs.modules[x] = 1

ii = leaves_list(hs.linkage)

mod_reordered = hs.modules.iloc[ii]

module_scores = hs.calculate_module_scores()
adata_result = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/PRPC_hs.adata_velocity_pseudotime.h5ad"
)

adata_result.obs["Module1"] = module_scores.loc[adata_result.obs.index, 1].values
adata_result.obs["Module2"] = module_scores.loc[adata_result.obs.index, 2].values
adata_result.obs["Module3"] = module_scores.loc[adata_result.obs.index, 3].values

adata_result.obs["Module"] = adata_result.obs[
    [
        "Module1",
        "Module2",
        "Module3",
    ]
].idxmax(axis=1)

module_list = [
    "Module1",
    "Module2",
    "Module3",
]
cols = px.colors.qualitative.D3

adata = adata[adata_result.obs.index,]
n = 20
adata_result.obs["velocity_pseudotime_cut"] = pd.cut(
    adata_result.obs.velocity_pseudotime, n
)


for module in [0, 1, 2]:
    plt.clf()
    z_score_global = [0] * n
    ax = plt.axes([0, 0, 1, 1], frameon=False)
    # Then we disable our xaxis and yaxis completely. If we just say plt.axis('off'),
    # they are still used in the computation of the image padding.
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # Even though our axes (plot region) are set to cover the whole image with [0,0,1,1],
    # by default they leave padding between the plotted data and the frame. We use tigher=True
    # to make sure the data gets scaled to the full extents of the axes.
    plt.autoscale(tight=True)

    for gene in mod_reordered[mod_reordered == module + 1].index[:20]:
        print(module)
        print(gene)
        y = adata[:, gene].X.toarray()
        y = [x[0] for x in y]
        n_convolve = 100
        weights = np.ones(n_convolve) / n_convolve
        y_hat = np.convolve(y, weights, mode="same")
        adata_result.obs["y_hat"] = y_hat
        df = (
            adata_result.obs.loc[
                :, ["velocity_pseudotime", "y_hat", "velocity_pseudotime_cut"]
            ]
            .groupby("velocity_pseudotime_cut")
            .mean()
        )
        y_pred = df.y_hat
        x = df.velocity_pseudotime
        z_scores = (y_pred - y_pred.mean()) / (y_pred.std())
        z_score_global = z_score_global + z_scores
        plt.plot(x, z_scores, color=cols[module], linewidth=0.8)
    z_score_global = z_score_global / len(
        mod_reordered[mod_reordered == module + 1].index[:20]
    )
    plt.plot(x, z_score_global, color="black", linewidth=2)
    plt.xticks([])
    plt.yticks([])
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.yticks(rotation=0)
    plt.savefig(
        output_file_path + "PRPC_gene_module_" + str(module + 1) + "_gene_trend.tiff",
        bbox_inches="tight",
        transparent=True,
        dpi=600,
    )
    plt.clf()


###
module=0
plt.clf()
z_score_global = [0] * n
ax = plt.axes([0, 0, 1, 1], frameon=False)
# Then we disable our xaxis and yaxis completely. If we just say plt.axis('off'),
# they are still used in the computation of the image padding.
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# Even though our axes (plot region) are set to cover the whole image with [0,0,1,1],
# by default they leave padding between the plotted data and the frame. We use tigher=True
# to make sure the data gets scaled to the full extents of the axes.
plt.autoscale(tight=True)

for gene in ["C21orf58", "MIS18BP1", "LINC01572", "SMC4", "RTKN2", "ASPM", "KIF14", "APOLD1", "MKI67", "TOP2A"] + ["LINC01414", "PCDH7", "NALCN-AS1", "KIAA1217", "ROBO2", "CNTN5", "MEG8", "FIGN", "PLEKHG1", "ESRRG"]:
    print(module)
    print(gene)
    y = adata[:, gene].X.toarray()
    y = [x[0] for x in y]
    n_convolve = 100
    weights = np.ones(n_convolve) / n_convolve
    y_hat = np.convolve(y, weights, mode="same")
    adata_result.obs["y_hat"] = y_hat
    df = (
        adata_result.obs.loc[
            :, ["velocity_pseudotime", "y_hat", "velocity_pseudotime_cut"]
        ]
        .groupby("velocity_pseudotime_cut")
        .mean()
    )
    y_pred = df.y_hat
    x = df.velocity_pseudotime
    z_scores = (y_pred - y_pred.mean()) / (y_pred.std())
    z_score_global = z_score_global + z_scores
    plt.plot(x, z_scores, color=cols[module], linewidth=0.8)
z_score_global = z_score_global / len(
    mod_reordered[mod_reordered == module + 1].index[:20]
)
plt.plot(x, z_score_global, color="black", linewidth=2)
plt.xticks([])
plt.yticks([])
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.yticks(rotation=0)
plt.savefig(
    output_file_path + "PRPC_gene_module_" + str(module + 1) + "_gene_trend.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
)
plt.clf()