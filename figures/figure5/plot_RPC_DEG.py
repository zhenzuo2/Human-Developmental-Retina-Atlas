import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({"font.size": 15})

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
adata_result = adata_result[adata_result.obs.majorclass.isin(["NRPC", "PRPC"])]
sc.pp.highly_variable_genes(
    adata_result, flavor="seurat_v3", n_top_genes=2000, subset=False
)
adata_result.var.highly_variable[
    [
        "CADPS2",
        "CHSY3",
        "COL25A1",
        "PDLIM5",
        "PTPRQ",
        "RMST",
        "EFNA5",
        "HMX1",
        "C14orf39",
        "ANOS1",
        "BAMBI",
        "CDH7",
        "HYDIN",
        "KCNIP4",
        "LINC01995",
        "OPCML",
        "SNCAIP",
        "PCDH9",
        "PTPRM",
        "ROBO2",
        "ANK3",
        "ZFPM2",
        "SLC22A23",
        "SLCO5A1",
        "AGBL4",
        "DCC",
        "L3MBTL4",
        "TTYH2",
        "WIF1",
        "EPHA3",
        "UNC5C",
        "WNT5B",
        "CALM1",
        "CGNL1",
        "SAT1",
        "TLL1",
        "CDH12",
        "ERBB4",
        "MAMDC2",
        "PPFIA2",
        "RARB",
        "CDH20",
        "SPOCK1",
        "MIR99AHG",
        "SLC8A1",
        "FBXL14",
        "P3H2",
        "SOX6",
        "TMTC1",
        "ARHGAP31",
        "OAF",
        "DIRC3",
        "TENM2",
    ]
] = True
adata_result = adata_result[:, adata_result.var.highly_variable]
sc.pp.normalize_total(adata_result)
sc.pp.scale(adata_result)
marker_genes_dict = {
    "Peripheral": [
        "CADPS2",
        "CHSY3",
        "COL25A1",
        "PDLIM5",
        "PTPRQ",
        "RMST",
        "EFNA5",
        "HMX1",
        "C14orf39",
        "ANOS1",
        "BAMBI",
        "CDH7",
        "HYDIN",
        "KCNIP4",
        "LINC01995",
        "OPCML",
        "SNCAIP",
        "PCDH9",
        "PTPRM",
        "ROBO2",
        "ANK3",
        "ZFPM2",
        "SLC22A23",
        "SLCO5A1",
        "AGBL4",
        "DCC",
        "L3MBTL4",
        "TTYH2",
        "WIF1",
    ],
    "Macula": [
        "EPHA3",
        "UNC5C",
        "WNT5B",
        "CALM1",
        "CGNL1",
        "SAT1",
        "TLL1",
        "CDH12",
        "ERBB4",
        "MAMDC2",
        "PPFIA2",
        "RARB",
        "CDH20",
        "SPOCK1",
        "MIR99AHG",
        "SLC8A1",
        "FBXL14",
        "P3H2",
        "SOX6",
        "TMTC1",
        "ARHGAP31",
        "OAF",
        "DIRC3",
        "TENM2",
    ],
}
adata_result.obs["Weeks"] = adata_result.obs.Days.map(
    {
        70: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW16",
        103: "PCW16",
        116: "PCW16",
        120: "PCW16",
        136: "PCW20",
        137: "PCW20",
        141: "PCW20",
        142: "PCW20",
        162: "PCW23",
        165: "PCW23",
    }
)
adata_result = adata_result[
    adata_result.obs.Weeks.isin(["PCW10", "PCW13", "PCW16", "PCW20"])
]


sc.pl.heatmap(
    adata_result,
    marker_genes_dict,
    groupby=[
        "Region",
        "Weeks",
    ],
    dendrogram=False,
    cmap="seismic",
    vmax=4,
    vmin=-4.5,
    show_gene_labels=True,
)
fig = plt.gcf()
fig.set_size_inches(10, 4)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/sc.pl.heatmap_RPC.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()
