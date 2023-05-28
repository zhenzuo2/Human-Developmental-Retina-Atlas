#!/usr/bin/env python
# coding: utf-8

import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sc.set_figure_params(scanpy=True, dpi=200)
colors = [
    "#440154FF",
    "#450558FF",
    "#46085CFF",
    "#470D60FF",
    "#471063FF",
    "#481467FF",
    "#481769FF",
    "#481B6DFF",
    "#481E70FF",
    "#482173FF",
    "#482576FF",
    "#482878FF",
    "#472C7AFF",
    "#472F7CFF",
    "#46327EFF",
    "#453581FF",
    "#453882FF",
    "#443B84FF",
    "#433E85FF",
    "#424186FF",
    "#404587FF",
    "#3F4788FF",
    "#3E4A89FF",
    "#3D4D8AFF",
    "#3C508BFF",
    "#3B528BFF",
    "#39558CFF",
    "#38598CFF",
    "#375B8DFF",
    "#355E8DFF",
    "#34608DFF",
    "#33638DFF",
    "#32658EFF",
    "#31688EFF",
    "#2F6B8EFF",
    "#2E6D8EFF",
    "#2D708EFF",
    "#2C718EFF",
    "#2B748EFF",
    "#2A768EFF",
    "#29798EFF",
    "#287C8EFF",
    "#277E8EFF",
    "#26818EFF",
    "#26828EFF",
    "#25858EFF",
    "#24878EFF",
    "#238A8DFF",
    "#228D8DFF",
    "#218F8DFF",
    "#20928CFF",
    "#20938CFF",
    "#1F968BFF",
    "#1F998AFF",
    "#1E9B8AFF",
    "#1F9E89FF",
    "#1FA088FF",
    "#1FA287FF",
    "#20A486FF",
    "#22A785FF",
    "#24AA83FF",
    "#25AC82FF",
    "#28AE80FF",
    "#2BB07FFF",
    "#2EB37CFF",
    "#31B67BFF",
    "#35B779FF",
    "#39BA76FF",
    "#3DBC74FF",
    "#41BE71FF",
    "#47C06FFF",
    "#4CC26CFF",
    "#51C56AFF",
    "#56C667FF",
    "#5BC863FF",
    "#61CA60FF",
    "#67CC5CFF",
    "#6DCD59FF",
    "#73D056FF",
    "#78D152FF",
    "#7FD34EFF",
    "#85D54AFF",
    "#8CD646FF",
    "#92D741FF",
    "#99D83DFF",
    "#A0DA39FF",
    "#A7DB35FF",
    "#ADDC30FF",
    "#B4DE2CFF",
    "#BBDE28FF",
    "#C2DF23FF",
    "#C9E020FF",
    "#D0E11CFF",
    "#D7E219FF",
    "#DDE318FF",
    "#E4E419FF",
    "#EBE51AFF",
    "#F1E51DFF",
    "#F7E620FF",
    # "#FDE725FF",
    "#5A5A5A",
]

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)

sc.pp.normalize_total(adata, exclude_highly_expressed=True)
sc.pp.log1p(adata)

sc.pl.umap(NRPC, color="leiden", legend_loc="on data")


# # AC
sc.pl.umap(NRPC, color="AC")

for clusters in set(NRPC.obs.leiden):
    print(clusters)
    clusters = [clusters]
    adata.obs["temp"] = np.nan
    adata.obs.loc[NRPC.obs[NRPC.obs.leiden.isin(clusters)].index, "temp"] = NRPC.obs[
        NRPC.obs.leiden.isin(clusters)
    ].AC
    adata.obs.loc[adata.obs.scpred_prediction == "AC", "temp"] = 1.01
    adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
    sc.pl.umap(
        adata,
        color="temp",
        palette=colors,
        legend_loc=None,
        frameon=False,
        title="",
    )

clusters = ["2", "13"]
print(clusters)
adata.obs["temp"] = np.nan
adata.obs.loc[NRPC.obs[NRPC.obs.leiden.isin(clusters)].index, "temp"] = NRPC.obs[
    NRPC.obs.leiden.isin(clusters)
].AC
adata.obs.loc[adata.obs.scpred_prediction == "AC", "temp"] = 1.01
adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
sc.pl.umap(
    adata, color="temp", palette=colors, legend_loc=None, frameon=False, title=""
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_leiden.svg",
    dpi=600,
    bbox_inches="tight",
)

sc.pl.umap(NRPC, color="NEUROD4", vmax="p95", vmin="p0")

sc.pl.umap(NRPC, color="TFAP2A", vmax="p95", vmin="p0")

sc.pl.umap(NRPC, color="TFAP2B", vmax="p95", vmin="p0")

sc.pl.umap(NRPC, color="ZEB2", vmax="p99", vmin="p0")

pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(clusters)].obs.index].obs,
        adata[adata.obs.scpred_prediction == "AC"].obs,
    ]
).to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_AC.csv"
)


# # HC
sc.pl.umap(NRPC, color="HC")

clusters = ["0", "1", "18"]
adata.obs["temp"] = np.nan
adata.obs.loc[NRPC.obs[NRPC.obs.leiden.isin(clusters)].index, "temp"] = NRPC.obs[
    NRPC.obs.leiden.isin(clusters)
].HC
adata.obs.loc[adata.obs.scpred_prediction.isin(["HC"]), "temp"] = 1.01
adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
sc.pl.umap(
    adata, color="temp", palette=colors, legend_loc=None, frameon=False, title=""
)

sc.pl.umap(NRPC, color="PROX1", vmax="p95", vmin="p0")

sc.pl.umap(NRPC, color="LHX1", vmax="p99", vmin="p0")

pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(clusters)].obs.index].obs,
        adata[adata.obs.scpred_prediction == "HC"].obs,
    ]
).to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_HC.csv"
)


# # BC
sc.pl.umap(NRPC, color="BC")

clusters = ["9", "10"]
adata.obs["temp"] = np.nan
adata.obs.loc[NRPC.obs[(NRPC.obs.leiden.isin(clusters))].index, "temp"] = NRPC.obs[
    NRPC.obs.leiden.isin(clusters)
].BC
adata.obs.loc[adata.obs.scpred_prediction == "BC", "temp"] = 1.01
adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
sc.pl.umap(
    adata, color="temp", palette=colors, legend_loc=None, frameon=False, title=""
)

sc.pl.umap(NRPC, color="OTX2", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="VSX2", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="ISL1", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="IRX6", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="GSG1", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="TRNP1", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="TMEM215", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="BHLHE22", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="FEZF2", vmax="p99", vmin="p0")

pd.concat(
    [
        adata[NRPC[(NRPC.obs.leiden.isin(clusters))].obs.index].obs,
        adata[adata.obs.scpred_prediction == "BC"].obs,
    ]
).to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_BC.csv"
)


# # Rod
sc.pl.umap(NRPC, color="Rod")

clusters = ["4"]
adata.obs["temp"] = np.nan
adata.obs.loc[NRPC.obs[(NRPC.obs.leiden.isin(clusters))].index, "temp"] = NRPC.obs[
    NRPC.obs.leiden.isin(clusters)
].Rod
adata.obs.loc[adata.obs.scpred_prediction == "Rod", "temp"] = 1.01
adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
sc.pl.umap(
    adata, color="temp", palette=colors, legend_loc=None, frameon=False, title=""
)

sc.pl.umap(NRPC, color="NR2E3", vmax="p99", vmin="p0")

pd.concat(
    [
        adata[NRPC[(NRPC.obs.leiden.isin(clusters))].obs.index].obs,
        adata[adata.obs.scpred_prediction == "Rod"].obs,
    ]
).to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_Rod.csv"
)


# # Cone
sc.pl.umap(NRPC, color="Cone")

clusters = ["17"]
adata.obs["temp"] = np.nan
adata.obs.loc[NRPC.obs[NRPC.obs.leiden.isin(clusters)].index, "temp"] = NRPC.obs[
    NRPC.obs.leiden.isin(clusters)
].Cone
adata.obs.loc[adata.obs.scpred_prediction == "Cone", "temp"] = 1.01
adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
sc.pl.umap(
    adata, color="temp", palette=colors, legend_loc=None, frameon=False, title=""
)

sc.pl.umap(NRPC, color="OTX2", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="CRX", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="PRDM1", vmax="p99", vmin="p0")

pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(clusters)].obs.index].obs,
        adata[adata.obs.scpred_prediction == "Cone"].obs,
    ]
).to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_Cone.csv"
)

# # RGC
sc.pl.umap(NRPC, color="RGC")

clusters = ["7", "14"]
adata.obs["temp"] = np.nan
adata.obs.loc[NRPC.obs[NRPC.obs.leiden.isin(clusters)].index, "temp"] = NRPC.obs[
    NRPC.obs.leiden.isin(clusters)
].RGC
adata.obs.loc[adata.obs.scpred_prediction == "RGC", "temp"] = 1.01
adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
sc.pl.umap(
    adata, color="temp", palette=colors, legend_loc=None, frameon=False, title=""
)

sc.pl.umap(NRPC, color="ATOH7", vmax="p99", vmin="p0")

sc.pl.umap(NRPC, color="POU4F2", vmax="p99", vmin="p0")

pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(clusters)].obs.index].obs,
        adata[adata.obs.scpred_prediction == "RGC"].obs,
    ]
).to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_RGC.csv"
)
