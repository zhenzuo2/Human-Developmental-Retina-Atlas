import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC.h5ad")
mv.velocity_graph(adata)
mv.latent_time(adata)

scv.pl.scatter(adata, color='NFIA', size=20,layer="ATAC",vmax = 1.7)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_ATAC.svg",
    dpi=600,
)

scv.pl.scatter(adata, color='NFIA', size=20,layer="Ms",vmax = 0.3)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_Ms.svg",
    dpi=600,
)

scv.pl.scatter(adata, color='NFIA', size=20,layer="Mu",vmax = 7)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_Mu.svg",
    dpi=600,
)

gene_list = ["NFIA"]
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
fig = plt.gcf()
fig.set_size_inches(10, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_state.svg",
    dpi=600,
)

gene_list = ["CBX1"]
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
fig = plt.gcf()
fig.set_size_inches(10, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/CBX1_state.svg",
    dpi=600,
)

gene_list = ["RPL10"]
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
fig = plt.gcf()
fig.set_size_inches(10, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/RPL10_state.svg",
    dpi=600,
)