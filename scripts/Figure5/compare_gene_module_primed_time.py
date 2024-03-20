import scanpy as sc
import multivelo as mv
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.legend_handler import HandlerTuple
import os
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import coo_matrix
from scipy.optimize import minimize
from scipy.spatial import KDTree
from sklearn.metrics import pairwise_distances
from sklearn.mixture import GaussianMixture
from scanpy import Neighbors
import scvelo as scv
import pandas as pd
import seaborn as sns
from numba import jit
from numba.typed import List
from tqdm.auto import tqdm
from joblib import Parallel, delayed


def my_dynamic_plot(
    adata,
    genes,
    by="expression",
    color_by="state",
    gene_time=True,
    axis_on=True,
    frame_on=True,
    show_anchors=True,
    show_switches=True,
    downsample=1,
    full_range=False,
    figsize=None,
    pointsize=2,
    linewidth=1.5,
    cmap="coolwarm",
):
    from pandas.api.types import is_numeric_dtype, is_categorical_dtype
    if by not in ["expression", "velocity"]:
        raise ValueError('"by" must be either "expression" or "velocity".')
    if by == "velocity":
        show_switches = False
    if color_by == "state":
        types = [0, 1, 2, 3]
        colors = ['#4169E1', '#008000', '#FF4500', '#800080']
    elif color_by in adata.obs and is_numeric_dtype(adata.obs[color_by]):
        types = None
        colors = adata.obs[color_by].values
    elif (
        color_by in adata.obs
        and is_categorical_dtype(adata.obs[color_by])
        and color_by + "_colors" in adata.uns.keys()
    ):
        types = adata.obs[color_by].cat.categories
        colors = adata.uns[f"{color_by}_colors"]
    else:
        raise ValueError(
            "Currently, color key must be a single string of "
            "either numerical or categorical available in adata "
            "obs, and the colors of categories can be found in "
            "adata uns."
        )
    downsample = np.clip(int(downsample), 1, 10)
    genes = np.array(genes)
    missing_genes = genes[~np.isin(genes, adata.var_names)]
    if len(missing_genes) > 0:
        print(f"{missing_genes} not found")
    genes = genes[np.isin(genes, adata.var_names)]
    gn = len(genes)
    if gn == 0:
        return
    if not gene_time:
        show_anchors = False
        latent_time = np.array(adata.obs["latent_time"])
        time_window = latent_time // 0.05
        time_window = time_window.astype(int)
        time_window[time_window == 20] = 19
    if "velo_s_params" in adata.uns.keys() and "outlier" in adata.uns["velo_s_params"]:
        outlier = adata.uns["velo_s_params"]["outlier"]
    else:
        outlier = 99
    fig, axs = plt.subplots(
        gn, 3, squeeze=False, figsize=(10, 2.3 * gn) if figsize is None else figsize
    )
    fig.patch.set_facecolor("white")
    for row, gene in enumerate(genes):
        u = adata[:, gene].layers["Mu" if by == "expression" else "velo_u"]
        s = adata[:, gene].layers["Ms" if by == "expression" else "velo_s"]
        c = adata[:, gene].layers["ATAC" if by == "expression" else "velo_chrom"]
        c = c.A if sparse.issparse(c) else c
        u = u.A if sparse.issparse(u) else u
        s = s.A if sparse.issparse(s) else s
        c, u, s = np.ravel(c), np.ravel(u), np.ravel(s)
        non_outlier = c <= np.percentile(c, outlier)
        non_outlier &= u <= np.percentile(u, outlier)
        non_outlier &= s <= np.percentile(s, outlier)
        c, u, s = c[non_outlier], u[non_outlier], s[non_outlier]
        time = np.array(adata[:, gene].layers["fit_t"] if gene_time else latent_time)
        if by == "velocity":
            time = np.reshape(time, (-1, 1))
            time = np.ravel(adata.obsp["_RNA_conn"].dot(time))
        time = time[non_outlier]
        if types is not None:
            for i in range(len(types)):
                if color_by == "state":
                    filt = adata[non_outlier, gene].layers["fit_state"] == types[i]
                else:
                    filt = adata[non_outlier, :].obs[color_by] == types[i]
                filt = np.ravel(filt)
                if np.sum(filt) > 0:
                    axs[row, 0].scatter(
                        time[filt][::downsample],
                        c[filt][::downsample],
                        s=pointsize,
                        c=colors[i],
                        alpha=0.6,
                    )
                    axs[row, 1].scatter(
                        time[filt][::downsample],
                        u[filt][::downsample],
                        s=pointsize,
                        c=colors[i],
                        alpha=0.6,
                    )
                    axs[row, 2].scatter(
                        time[filt][::downsample],
                        s[filt][::downsample],
                        s=pointsize,
                        c=colors[i],
                        alpha=0.6,
                    )
        else:
            axs[row, 0].scatter(
                time[::downsample],
                c[::downsample],
                s=pointsize,
                c=colors[non_outlier][::downsample],
                alpha=0.6,
                cmap=cmap,
            )
            axs[row, 1].scatter(
                time[::downsample],
                u[::downsample],
                s=pointsize,
                c=colors[non_outlier][::downsample],
                alpha=0.6,
                cmap=cmap,
            )
            axs[row, 2].scatter(
                time[::downsample],
                s[::downsample],
                s=pointsize,
                c=colors[non_outlier][::downsample],
                alpha=0.6,
                cmap=cmap,
            )
        if not gene_time:
            window_count = np.zeros(20)
            window_mean_c = np.zeros(20)
            window_mean_u = np.zeros(20)
            window_mean_s = np.zeros(20)
            for i in np.unique(time_window[non_outlier]):
                idx = time_window[non_outlier] == i
                window_count[i] = np.sum(idx)
                window_mean_c[i] = np.mean(c[idx])
                window_mean_u[i] = np.mean(u[idx])
                window_mean_s[i] = np.mean(s[idx])
            window_idx = np.where(window_count > 20)[0]
            axs[row, 0].plot(
                window_idx * 0.05 + 0.025,
                window_mean_c[window_idx],
                linewidth=linewidth,
                color="black",
                alpha=0.5,
            )
            axs[row, 1].plot(
                window_idx * 0.05 + 0.025,
                window_mean_u[window_idx],
                linewidth=linewidth,
                color="black",
                alpha=0.5,
            )
            axs[row, 2].plot(
                window_idx * 0.05 + 0.025,
                window_mean_s[window_idx],
                linewidth=linewidth,
                color="black",
                alpha=0.5,
            )
        if show_anchors:
            n_anchors = adata.uns["velo_s_params"]["t"]
            t_sw_array = np.array(
                [
                    adata[:, gene].var["fit_t_sw1"],
                    adata[:, gene].var["fit_t_sw2"],
                    adata[:, gene].var["fit_t_sw3"],
                ]
            )
            t_sw_array = t_sw_array[t_sw_array < 20]
            min_idx = int(adata[:, gene].var["fit_anchor_min_idx"])
            max_idx = int(adata[:, gene].var["fit_anchor_max_idx"])
            old_t = np.linspace(0, 20, n_anchors)[min_idx : max_idx + 1]
            new_t = old_t - np.min(old_t)
            new_t = new_t * 20 / np.max(new_t)
            if by == "velocity" and not full_range:
                anchor_interval = 20 / (max_idx + 1 - min_idx)
                min_idx = int(adata[:, gene].var["fit_anchor_velo_min_idx"])
                max_idx = int(adata[:, gene].var["fit_anchor_velo_max_idx"])
                start = (
                    0
                    + (min_idx - adata[:, gene].var["fit_anchor_min_idx"])
                    * anchor_interval
                )
                end = (
                    20
                    + (max_idx - adata[:, gene].var["fit_anchor_max_idx"])
                    * anchor_interval
                )
                new_t = np.linspace(start, end, max_idx + 1 - min_idx)
            ax = axs[row, 0]
            a_c = (
                adata[:, gene]
                .varm["fit_anchor_c" if by == "expression" else "fit_anchor_c_velo"]
                .ravel()[min_idx : max_idx + 1]
            )
            if show_switches:
                for t_sw in t_sw_array:
                    if t_sw > 0:
                        ax.vlines(
                            t_sw,
                            np.min(c),
                            np.max(c),
                            colors="black",
                            linestyles="dashed",
                            alpha=0.5,
                        )
            ax.plot(
                new_t[0 : new_t.shape[0]],
                a_c,
                linewidth=linewidth,
                color="black",
                alpha=0.5,
            )
            ax = axs[row, 1]
            a_u = (
                adata[:, gene]
                .varm["fit_anchor_u" if by == "expression" else "fit_anchor_u_velo"]
                .ravel()[min_idx : max_idx + 1]
            )
            if show_switches:
                for t_sw in t_sw_array:
                    if t_sw > 0:
                        ax.vlines(
                            t_sw,
                            np.min(u),
                            np.max(u),
                            colors="black",
                            linestyles="dashed",
                            alpha=0.5,
                        )
            ax.plot(
                new_t[0 : new_t.shape[0]],
                a_u,
                linewidth=linewidth,
                color="black",
                alpha=0.5,
            )
            ax = axs[row, 2]
            a_s = (
                adata[:, gene]
                .varm["fit_anchor_s" if by == "expression" else "fit_anchor_s_velo"]
                .ravel()[min_idx : max_idx + 1]
            )
            if show_switches:
                for t_sw in t_sw_array:
                    if t_sw > 0:
                        ax.vlines(
                            t_sw,
                            np.min(s),
                            np.max(s),
                            colors="black",
                            linestyles="dashed",
                            alpha=0.5,
                        )
            ax.plot(
                new_t[0 : new_t.shape[0]],
                a_s,
                linewidth=linewidth,
                color="black",
                alpha=0.5,
            )
        axs[row, 0].set_title(
            f"{gene} ATAC" if by == "expression" else f"{gene} chromatin velocity"
        )
        axs[row, 0].set_xlabel("t" if by == "expression" else "~t")
        axs[row, 0].set_ylabel("c" if by == "expression" else "dc/dt")
        axs[row, 1].set_title(
            f"{gene} unspliced" + ("" if by == "expression" else " velocity")
        )
        axs[row, 1].set_xlabel("t" if by == "expression" else "~t")
        axs[row, 1].set_ylabel("u" if by == "expression" else "du/dt")
        axs[row, 2].set_title(
            f"{gene} spliced" + ("" if by == "expression" else " velocity")
        )
        axs[row, 2].set_xlabel("t" if by == "expression" else "~t")
        axs[row, 2].set_ylabel("s" if by == "expression" else "ds/dt")
        for j in range(3):
            ax = axs[row, j]
            if not axis_on:
                ax.xaxis.set_ticks_position("none")
                ax.yaxis.set_ticks_position("none")
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
            if not frame_on:
                ax.xaxis.set_ticks_position("none")
                ax.yaxis.set_ticks_position("none")
                ax.set_frame_on(False)
    fig.tight_layout()


adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/PRPC_full.h5ad"
)
adata.obs.to_csv("/storage/chentemp/zz4/adult_dev_compare/temp/PRPC_full.obs.csv")
adata.var.to_csv("/storage/chentemp/zz4/adult_dev_compare/temp/PRPC_full.var.csv")
#adata = adata[adata.obs.majorclass.isin(["PRPC"])]

plt.clf()
my_dynamic_plot(adata, ["ABCA1"])
fig = plt.gcf()
fig.set_size_inches(12, 3)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/ABCA1_PRPC.my_dynamic_plot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

plt.clf()
my_dynamic_plot(adata, ["ANP32B"])
fig = plt.gcf()
fig.set_size_inches(12, 3)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/ANP32B_PRPC.my_dynamic_plot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

plt.clf()
my_dynamic_plot(adata, ["CNTN1"])
fig = plt.gcf()
fig.set_size_inches(12, 3)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/CNTN1_PRPC.my_dynamic_plot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

modules = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/mod_reordered.csv"
)
modules.index = modules["Unnamed: 0"].values

def get_summary_coupled_on(adata):
    t_sw = (
        adata.var[["fit_t_sw1", "fit_t_sw2", "fit_t_sw3"]].copy()
    )
    t_sw = t_sw.mask(t_sw > 20, 20)
    t_sw = t_sw.mask(t_sw < 0)
    t_sw["interval 1"] = t_sw["fit_t_sw1"]
    t_sw["t_sw2 - t_sw1"] = t_sw["fit_t_sw2"] - t_sw["fit_t_sw1"]
    t_sw["t_sw3 - t_sw2"] = t_sw["fit_t_sw3"] - t_sw["fit_t_sw2"]
    t_sw["20 - t_sw3"] = 20 - t_sw["fit_t_sw3"]
    t_sw = t_sw.mask(t_sw <= 0)
    t_sw = t_sw.mask(t_sw > 20)
    t_sw.columns = pd.Index(
        [
            "time 1",
            "time 2",
            "time 3",
            "primed",
            "coupled-on",
            "decoupled",
            "coupled-off",
        ]
    )
    t_sw = t_sw[["primed", "coupled-on", "decoupled", "coupled-off"]]
    t_sw = t_sw / 20
    t_sw.to_csv("/storage/chentemp/zz4/adult_dev_compare/temp/t_sw.csv")
    a = t_sw["coupled-on"][~np.isnan(t_sw["coupled-on"])]
    return a

def get_summary_decoupled(adata):
    t_sw = (
        adata.var[["fit_t_sw1", "fit_t_sw2", "fit_t_sw3"]].copy()
    )
    t_sw = t_sw.mask(t_sw > 20, 20)
    t_sw = t_sw.mask(t_sw < 0)
    t_sw["interval 1"] = t_sw["fit_t_sw1"]
    t_sw["t_sw2 - t_sw1"] = t_sw["fit_t_sw2"] - t_sw["fit_t_sw1"]
    t_sw["t_sw3 - t_sw2"] = t_sw["fit_t_sw3"] - t_sw["fit_t_sw2"]
    t_sw["20 - t_sw3"] = 20 - t_sw["fit_t_sw3"]
    t_sw = t_sw.mask(t_sw <= 0)
    t_sw = t_sw.mask(t_sw > 20)
    t_sw.columns = pd.Index(
        [
            "time 1",
            "time 2",
            "time 3",
            "primed",
            "coupled-on",
            "decoupled",
            "coupled-off",
        ]
    )
    t_sw = t_sw[["primed", "coupled-on", "decoupled", "coupled-off"]]
    t_sw = t_sw / 20
    a = t_sw["decoupled"][~np.isnan(t_sw["decoupled"])]
    return a

df = pd.DataFrame(adata.var)
df["Module"] = ""
for gene in adata.var.index:
    if gene in modules.index:
        df.loc[gene, "Module"] = modules.loc[gene, "Module"]
    else:
        df.loc[gene, "Module"] = 0

###
res = get_summary_coupled_on(adata)
df = df.loc[res.index, :]
for gene in df.index:
    if gene in res.index:
        df.loc[gene, "coupled_on"] = res[gene]

df["Module"] = df["Module"].map(
    {
        1: "Module 1: Early differentiation",
        #2: "Module 2: In Transition",
        3: "Module 3: Late differentiation",
    }
)
df["Module"] = df["Module"].astype("category")
df["coupled_on"] = df["coupled_on"] * 100
df.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Compare Gene Coupled-on Time Percentage.csv")
plt.rcParams.update({"font.size": 15})
plt.clf()
ax = sns.boxplot(
    data=df,
    x="Module",
    y="coupled_on",
    boxprops={"alpha": 0.7},
    showfliers=False,
    palette=["#1f77b4", "#2ca02c"],
)
sns.stripplot(
    data=df,
    x="Module",
    y="coupled_on",
    hue="Module",
    palette=["#1f77b4", "#2ca02c"],
)
handles, labels = ax.get_legend_handles_labels()
plt.legend(bbox_to_anchor=(1.02, 0.21), loc="upper left", borderaxespad=0)
plt.xticks(rotation=45)  # Rotates X-Axis Ticks by 45-degrees
plt.ylabel("Gene Coupled-on Time Percentage")
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
fig = plt.gcf()
fig.set_size_inches(8, 8)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Compare Gene Coupled-on Time Percentage.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
###
df = pd.DataFrame(adata.var)
df["Module"] = ""
for gene in adata.var.index:
    if gene in modules.index:
        df.loc[gene, "Module"] = modules.loc[gene, "Module"]
    else:
        df.loc[gene, "Module"] = 0
res = get_summary_decoupled(adata)
df = df.loc[res.index, :]
for gene in df.index:
    if gene in res.index:
        df.loc[gene, "decoupled"] = res[gene]

df["Module"] = df["Module"].map(
    {
        1: "Module 1: Early differentiation",
        #2: "Module 2: In Transition",
        3: "Module 3: Late differentiation",
    }
)
df["Module"] = df["Module"].astype("category")
df["decoupled"] = df["decoupled"] * 100
df.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Compare Gene Decoupled Time Percentage.csv")
plt.rcParams.update({"font.size": 15})
plt.clf()
ax = sns.boxplot(
    data=df,
    x="Module",
    y="decoupled",
    boxprops={"alpha": 0.7},
    showfliers=False,
    palette=["#1f77b4", "#2ca02c"],
)
sns.stripplot(
    data=df,
    x="Module",
    y="decoupled",
    hue="Module",
    palette=["#1f77b4", "#2ca02c"],
)
handles, labels = ax.get_legend_handles_labels()
plt.legend(bbox_to_anchor=(1.02, 0.21), loc="upper left", borderaxespad=0)
plt.xticks(rotation=45)  # Rotates X-Axis Ticks by 45-degrees
plt.ylabel("Gene Decoupled Time Percentage")
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
fig = plt.gcf()
fig.set_size_inches(8, 8)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Compare Gene Decoupled Time Percentage.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
###
plt.clf()
my_dynamic_plot(adata, ["NDC80", "DEPDC1B", "SKA1", "AURKB", "HMMR", "KIF11",])
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Module1_PRPC.my_dynamic_plot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

plt.clf()
my_dynamic_plot(adata, ["ABCG1", "ATP2B2", "SAT1", "SEMA5B", "ST3GAL1", "WIF1"])
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Module3_PRPC.my_dynamic_plot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

plt.clf()
my_dynamic_plot(adata, ["RPL10", "RPL12", "RPL15", "RPL17", 
"RPL18", "RPL21"])
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Module2_RPL_PRPC.my_dynamic_plot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

plt.clf()
my_dynamic_plot(adata, ["RPS14", "RPS15A", "RPS17", "RPS20", "RPS27", "RPS29"])
fig = plt.gcf()
fig.set_size_inches(12, 8)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/Module2_RPS_PRPC.my_dynamic_plot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

def get_summary(adata):
    t_sw = (
        adata#[:, adata.var["velo_s_genes"]]
        .var[["fit_t_sw1", "fit_t_sw2", "fit_t_sw3"]]
        .copy()
    )
    t_sw = t_sw.mask(t_sw > 20, 20)
    t_sw = t_sw.mask(t_sw < 0)
    t_sw["interval 1"] = t_sw["fit_t_sw1"]
    t_sw["t_sw2 - t_sw1"] = t_sw["fit_t_sw2"] - t_sw["fit_t_sw1"]
    t_sw["t_sw3 - t_sw2"] = t_sw["fit_t_sw3"] - t_sw["fit_t_sw2"]
    t_sw["20 - t_sw3"] = 20 - t_sw["fit_t_sw3"]
    t_sw = t_sw.mask(t_sw <= 0)
    t_sw = t_sw.mask(t_sw > 20)
    t_sw.columns = pd.Index(
        [
            "time 1",
            "time 2",
            "time 3",
            "primed",
            "coupled-on",
            "decoupled",
            "coupled-off",
        ]
    )
    t_sw = t_sw[["primed", "coupled-on", "decoupled", "coupled-off"]]
    t_sw = t_sw / 20
    a = t_sw["primed"][~np.isnan(t_sw["primed"])]
    b = t_sw["coupled-on"][~np.isnan(t_sw["coupled-on"])]
    c = t_sw["decoupled"][~np.isnan(t_sw["decoupled"])]
    d = t_sw["coupled-off"][~np.isnan(t_sw["coupled-off"])]
    return (
        list(a.values) + list(b.values) + list(c.values) + list(d.values),
        len(a) * ["Primed"]
        + len(b) * ["Coupled-on"]
        + len(c) * ["Decoupled"]
        + len(d) * ["Coupled-off"],
    )


a, b = get_summary(adata)
d = {"col1": [x * 100 for x in a], "col2": b}
pd.DataFrame(data=d).to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/switch_time_summary.csv")
# Set up colors for the stacked bars
colors = ['#4169E1', '#008000', '#FF4500', '#800080']
plt.clf()
ax = sns.boxplot(d, x="col2", y="col1", palette=colors, showfliers=False)

# Add labels and title
plt.xlabel("")
plt.ylabel("")
plt.title("")
plt.legend("", frameon=False)
ax.yaxis.set_major_formatter(mtick.PercentFormatter())

fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks(rotation=90)  # Rotates X-Axis Ticks by 45-degrees
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/switch_time_summary.png",
    dpi=600,
    bbox_inches="tight",
)
