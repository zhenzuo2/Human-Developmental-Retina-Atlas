import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy.spatial import ConvexHull
from sklearn.neighbors import NearestNeighbors

matplotlib.rcParams.update({"font.size": 22})

PRPC = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/PRPC_full.h5ad"
)
NRPC = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/NRPC.h5ad"
)
AC = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/AC.h5ad"
)
BC = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/BC.h5ad"
)
Rod = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/Rod.h5ad"
)
Cone = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/Cone.h5ad"
)
RGC = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/RGC.h5ad"
)
HC = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/HC.h5ad"
)

def get_summary(adata):
    t_sw = (
        adata[:, adata.var["velo_s_genes"]]
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
    a = np.median(a)
    b = t_sw["coupled-on"][~np.isnan(t_sw["coupled-on"])]
    b = np.median(b)
    c = t_sw["decoupled"][~np.isnan(t_sw["decoupled"])]
    c = np.median(c)
    d = t_sw["coupled-off"][~np.isnan(t_sw["coupled-off"])]
    d = np.median(d)
    return [a, b, c, d]


# Generate random data for three time points
time_points = [
    "PRPC",
    "NRPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
]
observations = [
    get_summary(PRPC),
    get_summary(NRPC),
    get_summary(RGC),
    get_summary(Cone),
    get_summary(HC),
    get_summary(AC),
    get_summary(Rod),
    get_summary(BC),
]
observations = observations / np.sum(observations, axis=1)[:, np.newaxis]
np.savetxt('/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/switch_time_summary.csv', observations.T, delimiter=',')

# Set up colors for the stacked bars
colors = ["#4169E1", "#008000", "#FF4500", "#800080"]

# Plot the stacked bars

bottom = np.zeros(len(time_points))
plt.clf()

for i, obs in enumerate(observations.T):
    plt.bar(time_points, obs, bottom=bottom, color=colors[i])
    bottom += obs

# Add labels and title
plt.xlabel("")
plt.ylabel("")
plt.title("")
# Show the legend
plt.legend(
    ["Primed", "Coupled-on", "Decoupled", "Coupled-off"],
    bbox_to_anchor=(1.04, 1),
    loc="upper left",
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks(
    ticks=[0, 1, 2, 3, 4, 5, 6, 7],
    labels=[
        "PRPC",
        "NRPC",
        "RGC",
        "Cone",
        "HC",
        "AC",
        "Rod",
        "BC",
    ],
    rotation=45,
)
from matplotlib.ticker import FuncFormatter
def percentage_formatter(x, pos):
    return f'{x * 100:.0f}%'
plt.gca().yaxis.set_major_formatter(FuncFormatter(percentage_formatter))
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/switch_time_summary.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


def get_summary(adata):
    genes = adata.var_names
    fit_model = (
        adata[
            :,
            (adata.var["fit_direction"] == "complete")
            & np.isin(adata.var_names, genes),
        ]
        .var["fit_model"]
        .values
    )
    fit_direction = adata[:, genes].var["fit_direction"].values
    data = [
        np.sum(fit_direction == "on"),
        np.sum(fit_direction == "off"),
        np.sum(fit_model == 1),
        np.sum(fit_model == 2),
    ]
    index = ["induction", "repression", "Model 1", "Model 2"]
    index = [x for i, x in enumerate(index) if data[i] > 0]
    data = [x for x in data if x > 0]
    df = pd.DataFrame({"data": data}, index=index)
    return df


observations = [
    get_summary(PRPC),
    get_summary(NRPC),
    get_summary(RGC),
    get_summary(Cone),
    get_summary(HC),
    get_summary(AC),
    get_summary(Rod),
    get_summary(BC),
]
res = np.zeros((8, 4))
for i in range(8):
    j = 0
    for x in ["induction", "repression", "Model 1", "Model 2"]:
        if x in observations[i].index:
            res[i, j] = observations[i].loc[x, "data"]
        j = j + 1


res = res / np.sum(res, axis=1)[:, np.newaxis]

# Set up colors for the stacked bars
colors = ["#4169E1", "#008000", "#FF4500", "#800080"]

# Plot the stacked bars

bottom = np.zeros(len(time_points))
plt.clf()
for i, obs in enumerate(res.T):
    plt.bar(time_points, obs, bottom=bottom, color=colors[i])
    bottom += obs

# Add labels and title
plt.xlabel("")
plt.ylabel("")
plt.title("")

# Show the legend
plt.legend(
    ["Induction", "Repression", "Model 1", "Model 2"],
    labelcolor=colors,
    bbox_to_anchor=(1.04, 1),
    loc="upper left",
)

fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks(
    ticks=[0, 1, 2, 3, 4, 5, 6, 7],
    labels=[
        "PRPC",
        "NRPC",
        "RGC",
        "Cone",
        "HC",
        "AC",
        "Rod",
        "BC",
    ],
    rotation=45,
)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/model_summary.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
