import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
NRPC = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC = NRPC[NRPC.obs.leiden!="12"]
NRPC.obs["subclass"] = np.nan

sc.pp.neighbors(NRPC, use_rep="X_scVI")
sc.tl.leiden(NRPC, resolution=3)

NRPC.obs['temp'] = NRPC.obs.leiden.map({'21': 'TBD', '27': 'OFF Cone BC','14':"RBC/On Cone BC"})

sc.pl.umap(
    NRPC,
    size=40,
    legend_loc=None,
    colorbar_loc=None,
    color=["temp"],
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/NRPC_BC_annotation.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'

legend_dict = {
    "BC Precursor": "#1f77b4",
    "OFF Cone BC (NRPC)": "#ff7f0e",
    "OFF Cone BC": "#2ca02c",
    "ON Cone BC": "#d62728",
    "RBC": "#9467bd",
    "RBC/On Cone BC (NRPC)": "#8c564b",
    "TBD (NRPC)": "#e377c2"
}

lines = [mlines.Line2D([], [], color=color, marker='o', markersize=25, 
            markerfacecolor=color, label=label, linestyle='None') for label, color in legend_dict.items()]

plt.figure(figsize=(10,5))
plt.legend(handles=lines, loc='center', ncol=2, fontsize=20)
plt.axis('off')
plt.savefig("legend.svg")