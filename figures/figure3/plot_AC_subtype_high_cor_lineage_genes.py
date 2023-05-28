import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
import joblib
import cellrank as cr

sc.set_figure_params(dpi_save=600, transparent=True)

output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/"
g = joblib.load(
    "/storage/singlecell/zz4/fetal_snakemake/results/cellrank_result/AC_g.h5ad",
)
adata = g.adata
adata.obs.loc[adata.obs.leiden == "12", "majorclass"] = "AC Precursor"
sc.pl.umap(adata, color="majorclass")

adata.obs["GABAergic"] = [
    x[0] for x in adata.obsm["to_terminal_states"]["GABAergic"].tolist()
]
adata.obs["temp"] = np.nan
adata.obs.loc[adata.obs.GABAergic > 0.83, "temp"] = adata.obs.loc[
    adata.obs.GABAergic > 0.83, "GABAergic"
]
sc.pl.umap(adata, color="temp", colorbar_loc=None, frameon=False, title="", size=50)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "GABAergic_p.svg",
    bbox_inches="tight",
)
drivers = g.compute_lineage_drivers(lineages="GABAergic", return_drivers=True)
((drivers.GABAergic_qval < 0.01) & (drivers.GABAergic_ci_low >= 0.1)).value_counts()
top_genes = drivers[(drivers.GABAergic_qval < 0.01)].index[:146]
ax = scv.pl.heatmap(
    adata[adata.obs.GABAergic > 0.83],
    var_names=top_genes,
    sortby="GABAergic",
    col_color="majorclass",
    yticklabels=True,
    show=False,
)
genes = [x.get_text() for x in ax.ax_heatmap.get_yticklabels()][:50]
scv.pl.heatmap(
    adata[adata.obs.GABAergic > 0.83],
    var_names=genes,
    sortby="GABAergic",
    col_color="majorclass",
    font_scale=0.5,
    yticklabels=True,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "GABAergic_cor_genes.png",
    bbox_inches="tight",
)

adata.obs["Glycinergic"] = [
    x[0] for x in adata.obsm["to_terminal_states"]["Glycinergic"].tolist()
]
adata.obs["temp"] = np.nan
adata.obs.loc[adata.obs.Glycinergic > 0.07, "temp"] = adata.obs.loc[
    adata.obs.Glycinergic > 0.07, "Glycinergic"
]
sc.pl.umap(adata, color="temp", colorbar_loc=None, frameon=False, title="", size=50)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "Glycinergic_p.svg",
    bbox_inches="tight",
)
drivers = g.compute_lineage_drivers(lineages="Glycinergic", return_drivers=True)
drivers.sort_values(by="Glycinergic_corr", ascending=False)
((drivers.Glycinergic_qval < 0.01) & (drivers.Glycinergic_ci_low >= 0.2)).value_counts()
top_genes = drivers[
    (drivers.Glycinergic_qval < 0.01) & (drivers.Glycinergic_ci_low >= 0.1)
].index[:160]
ax = scv.pl.heatmap(
    adata[adata.obs.Glycinergic > 0.07],
    var_names=top_genes,
    sortby="Glycinergic",
    col_color="majorclass",
    yticklabels=True,
    show=False,
)
genes = [x.get_text() for x in ax.ax_heatmap.get_yticklabels()][15:50]
scv.pl.heatmap(
    adata[adata.obs.Glycinergic > 0.07],
    var_names=genes,
    sortby="Glycinergic",
    col_color="majorclass",
    font_scale=0.5,
    yticklabels=True,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "Glycinergic_cor_genes.png",
    bbox_inches="tight",
)

adata.obs["SACs"] = [x[0] for x in adata.obsm["to_terminal_states"]["SACs"].tolist()]
adata.obs["temp"] = np.nan
adata.obs.loc[adata.obs.SACs > 0.1, "temp"] = adata.obs.loc[
    adata.obs.SACs > 0.1, "SACs"
]
sc.pl.umap(adata, color="temp", colorbar_loc=None, frameon=False, title="", size=50)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "SACs_p.svg",
    bbox_inches="tight",
)
drivers = g.compute_lineage_drivers(lineages="SACs", return_drivers=True)
drivers.sort_values(by="SACs_corr", ascending=False)
((drivers.SACs_qval < 0.01) & (drivers.SACs_ci_low >= 0.1)).value_counts()
top_genes = drivers[(drivers.SACs_qval < 0.01) & (drivers.SACs_ci_low >= 0.1)].index[
    :150
]
ax = scv.pl.heatmap(
    adata[adata.obs.SACs > 0.02],
    var_names=top_genes,
    sortby="SACs",
    col_color="majorclass",
    yticklabels=True,
    show=False,
)
genes = [x.get_text() for x in ax.ax_heatmap.get_yticklabels()][3:50]
scv.pl.heatmap(
    adata[adata.obs.SACs > 0.02],
    var_names=genes,
    sortby="SACs",
    col_color="majorclass",
    font_scale=0.5,
    yticklabels=True,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "SACs_cor_genes.png",
    bbox_inches="tight",
)

adata.obs["dual ACs"] = [
    x[0] for x in adata.obsm["to_terminal_states"]["dual ACs"].tolist()
]
adata.obs["temp"] = np.nan
adata.obs.loc[adata.obs["dual ACs"] > 0.11, "temp"] = adata.obs.loc[
    adata.obs["dual ACs"] > 0.11, "dual ACs"
]
sc.pl.umap(adata, color="temp", colorbar_loc=None, frameon=False, title="", size=50)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "dual ACs_p.svg",
    bbox_inches="tight",
)
drivers = g.compute_lineage_drivers(lineages="dual ACs", return_drivers=True)
drivers.sort_values(by="dual ACs_corr", ascending=False)
((drivers["dual ACs_qval"] < 0.01) & (drivers["dual ACs_ci_low"] >= 0.1)).value_counts()
top_genes = drivers[
    (drivers["dual ACs_qval"] < 0.01) & (drivers["dual ACs_ci_low"] >= 0.1)
].index[:150]
ax = scv.pl.heatmap(
    adata[adata.obs["dual ACs"] > 0.11],
    var_names=top_genes,
    sortby="dual ACs",
    col_color="majorclass",
    yticklabels=True,
    show=False,
)
genes = [x.get_text() for x in ax.ax_heatmap.get_yticklabels()][3:50]
scv.pl.heatmap(
    adata[adata.obs["dual ACs"] > 0.11],
    var_names=genes,
    sortby="dual ACs",
    col_color="majorclass",
    font_scale=0.5,
    yticklabels=True,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "dual ACs_cor_genes.png",
    bbox_inches="tight",
)
