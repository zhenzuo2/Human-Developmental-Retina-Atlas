import scvelo as scv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import pandas as pd
import scanpy as sc
import seaborn as sns

matplotlib.rcParams.update({"font.size": 30})
adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
adata.obs["Weeks"] = adata.obs.Days.map(
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
adata.obs["Weeks"] = (
    adata.obs["Weeks"]
    .astype("category")
    .cat.reorder_categories(["PCW10", "PCW13", "PCW16", "PCW20", "PCW23"])
)
adata[adata.obs.majorclass.isin(["Rod", "Cone"])].obs.Weeks.value_counts()
genes = [
    "OTX2",
    "RCVRN",
    "AIPL1",
    "NRL",
    "CRX",
    "PDE6B",
    "NR2E3",
    "ROM1",
    "GNGT2",
    "GNAT1",
    "PDE6H",
]
for gene in genes:
    adata.obs[gene] = adata[:, gene].X.toarray()
    plt.clf()
    sns.violinplot(
        data=adata[adata.obs.majorclass.isin(["Rod", "Cone"])].obs,
        x="Weeks",
        y=gene,
        hue="Region",
        split=True,
        inner=None,
        palette={"Macula": "#F8766D", "Peripheral": "#00BFC4"},
    )
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    #plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.legend('',frameon=False)
    plt.xticks(rotation='vertical')
    plt.xlabel('')
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/sns.violinplot_Rod_Cone_"
        + gene
        + ".svg",
        dpi=600,
        bbox_inches="tight",
        transparent=True,
        backend="cairo",
    )
    plt.clf()

genes = ["GRM6", "LHX4", "PRDM8", "CABP5", "VSX1"]
for gene in genes:
    adata.obs[gene] = adata[:, gene].X.toarray()
    plt.clf()
    sns.violinplot(
        data=adata[adata.obs.majorclass.isin(["BC"])].obs,
        x="Weeks",
        y=gene,
        hue="Region",
        split=True,
        inner=None,
        palette={"Macula": "#F8766D", "Peripheral": "#00BFC4"},
    )
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.xticks(rotation='vertical')
    plt.xlabel('')
    #plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.legend('',frameon=False)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/sns.violinplot_BC_"
        + gene
        + ".svg",
        dpi=600,
        bbox_inches="tight",
        transparent=True,
        backend="cairo",
    )
    plt.clf()
