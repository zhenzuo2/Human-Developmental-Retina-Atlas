import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

pd.options.display.max_columns = None
scv.set_figure_params(dpi=600, dpi_save=600)

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.h5ad"
)
sc.pp.neighbors(adata_result, use_rep="X_scANVI")

TFs = [
    "ARID3A",
    "ARID5B",
    "ARNTL2",
    "ASCL1",
    "ATOH7",
    "BACH2",
    "BCL11A",
    "CCDC88A",
    "CDC5L",
    "CPEB1",
    "CREB5",
    "CTCF",
    "CUX2",
    "DACH1",
    "DPF1",
    "E2F1",
    "E2F2",
    "E2F3",
    "E2F7",
    "E2F8",
    "EGR1",
    "ELF2",
    "EPAS1",
    "ESR2",
    "ESRRG",
    "FEZF2",
    "FOS",
    "FOXM1",
    "FOXN3",
    "FOXN4",
    "FOXO1",
    "FOXP1",
    "FOXP2",
    "FOXP4",
    "GABPB1",
    "GLI2",
    "GLI3",
    "GLIS3",
    "HES1",
    "HES4",
    "HES5",
    "HES6",
    "HEY2",
    "HIF1A",
    "HIF3A",
    "HIVEP2",
    "HIVEP3",
    "HMGA1",
    "HMGA2",
    "HMGB1",
    "HMGB2",
    "HMGB3",
    "HMGXB4",
    "HMX1",
    "HOPX",
    "ID2",
    "ID3",
    "IKZF2",
    "IRX5",
    "IRX6",
    "JAZF1",
    "JUNB",
    "KLF13",
    "KLF6",
    "LCORL",
    "LEF1",
    "LHX9",
    "LIN28B",
    "MAEL",
    "MAZ",
    "MECOM",
    "MEF2C",
    "MEIS1",
    "MITF",
    "MLLT10",
    "MXD3",
    "MXI1",
    "MYB",
    "MYBL1",
    "MYBL2",
    "NCOA1",
    "NEUROD4",
    "NFATC2",
    "NFATC3",
    "NFIA",
    "NFIB",
    "NFIX",
    "NPAS2",
    "NPAS3",
    "NR2C2",
    "NR2F1",
    "NR2F2",
    "NR3C2",
    "NR4A1",
    "NR4A3",
    "NRF1",
    "NSD2",
    "ONECUT1",
    "OTX2",
    "OVOL2",
    "PBX1",
    "PBX3",
    "PBX4",
    "PKNOX2",
    "PLAG1",
    "PLAGL1",
    "POU2F1",
    "POU3F1",
    "PRDM5",
    "PROX1",
    "RARB",
    "RAX",
    "RFX2",
    "RORB",
    "SMAD2",
    "SMAD6",
    "SOX11",
    "SOX13",
    "SOX2",
    "SOX4",
    "SOX5",
    "SOX6",
    "SOX8",
    "SP3",
    "SREBF2",
    "ST18",
    "STAT3",
    "TBX20",
    "TBX5",
    "TCF12",
    "TCF4",
    "TCF7L1",
    "TCF7L2",
    "TEAD1",
    "TFDP1",
    "TFDP2",
    "THRB",
    "TOX",
    "TOX3",
    "TP73",
    "TRPS1",
    "TSHZ2",
    "TSHZ3",
    "YBX1",
    "ZBTB16",
    "ZBTB18",
    "ZBTB20",
    "ZBTB7C",
    "ZEB1",
    "ZFHX3",
    "ZFPM2",
    "ZIC1",
    "ZIC5",
    "ZNF100",
    "ZNF107",
    "ZNF257",
    "ZNF276",
    "ZNF331",
    "ZNF367",
    "ZNF385D",
    "ZNF423",
    "ZNF483",
    "ZNF492",
    "ZNF516",
    "ZNF519",
    "ZNF521",
    "ZNF536",
    "ZNF569",
    "ZNF66",
    "ZNF680",
    "ZNF682",
    "ZNF695",
    "ZNF704",
    "ZNF714",
    "ZNF718",
    "ZNF724",
    "ZNF726",
    "ZNF730",
    "ZNF732",
    "ZNF827",
    "ZNF83",
    "ZNF85",
    "ZNF850",
    "ZNF90",
    "ZNF93",
    "TADA2A",
    "ZHX3",
    "DACH2",
    "EOMES",
    "FIGLA",
    "FOXJ1",
    "HLF",
    "JUN",
    "KLF9",
    "MYOG",
    "NEUROD1",
    "NFAT5",
    "NFE2",
    "OTP",
    "PAX2",
    "PPARG",
    "PRDM8",
    "SMAD3",
    "SSRP1",
    "TAL2",
    "TEAD4",
    "TFAP2D",
    "TFEC",
    "VAX2",
    "ZNF438",
    "ZNF467",
    "NR3C1",
    "RORA",
    "FOSL1",
    "LYL1",
    "NRL",
    "RBPJ",
    "SALL3",
    "SALL4",
    "WT1",
    "ZFP91",
    "ZNF486",
]

adata_result = adata_result[:, [x for x in TFs if x in adata_result.var.index]]
adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
scv.pl.umap(adata_result, color="Days")

adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
scv.pl.umap(adata_result[adata_result.obs.Region == "Macula"], color="Days")

adata_result.obs["Days"] = adata_result.obs["Days"].astype(float)
scv.pl.umap(adata_result[adata_result.obs.Region == "Peripheral"], color="Days")

adata_result_Macula = adata_result[adata_result.obs.Region =="Macula"]
adata_result_Peripheral = adata_result[adata_result.obs.Region =="Peripheral"]
sc.pp.neighbors(adata_result_Macula, use_rep="X_scANVI")
sc.pp.neighbors(adata_result_Peripheral, use_rep="X_scANVI")

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

mv.velocity_graph(adata_result_Macula)
mv.latent_time(adata_result_Macula)

mv.velocity_graph(adata_result_Peripheral)
mv.latent_time(adata_result_Peripheral)

mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="Days",
    color_map="viridis",
    title="",
    save="velocity_embedding_stream.svg",
)

mv.velocity_embedding_stream(
    adata_result[adata_result.obs.Region =="Macula"],
    basis="umap",
    color="Days",
    color_map="viridis",
    title="",
    save="velocity_embedding_stream_Macula.svg",
)

mv.velocity_embedding_stream(
    adata_result[adata_result.obs.Region =="Peripheral"],
    basis="umap",
    color="Days",
    color_map="viridis",
    title="",
    save="velocity_embedding_stream_Peripheral.svg",
)


adata_result.uns["velo_s_graph"] = adata_result.uns["velo_s_norm_graph"]

scv.tl.velocity_confidence(adata_result, vkey="velo_s")
keys = "velo_s_length", "velo_s_confidence"
scv.pl.scatter(adata_result, c=keys, cmap="coolwarm", perc=[5, 95])

scv.tl.velocity_graph(adata_result)
scv.tl.velocity_confidence(adata_result)

keys = "velocity_length", "velocity_confidence"
scv.pl.scatter(adata_result, c=keys, cmap="coolwarm", perc=[5, 95])

keys = "velocity_length"
scv.pl.scatter(
    adata_result, c=keys, cmap="coolwarm", perc=[5, 95], save="velocity_length.svg"
)


def map_values(map_from, map_to, values):
    mapped_values = []
    for value in values:
        index = map_from.index(value)
        mapped_values.append(map_to[index])
    return mapped_values


adata_result.obs["Days_group"] = map_values(
    map_from=[70, 79, 87, 91, 100, 103, 116, 136, 137, 141, 142, 162, 165],
    map_to=[
        "FW10",
        "FW10",
        "FW13",
        "FW13",
        "FW13",
        "FW16",
        "FW16",
        "FW19",
        "FW19",
        "FW19",
        "FW19",
        "FW23",
        "FW23",
    ],
    values=adata_result.obs.Days,
)

import seaborn as sns
import pandas as pd

# create example DataFrame

# create grouped boxplot using Seaborn
ax = sns.boxplot(x="Days_group", y="latent_time", hue="Region", data=adata_result.obs)
ax.set_xlabel("Fetal Weeks")
ax.set_ylabel("Latent Time")
plt.savefig("boxplot.svg")

adata_result.obs.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv"
)
adata_result.var.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.var.csv"
)