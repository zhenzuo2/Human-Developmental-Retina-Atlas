import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/HC_subtype_NRPC.h5ad")
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/"

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
cells = []
for sb in set(adata.obs.subclass):
    cells = cells + [adata[adata.obs.subclass==sb].obs.index[0]]

for time in set(adata.obs["Weeks"]):
    for region in set(adata.obs.Region):
        adata.obs["temp"] = (
            (adata.obs.Weeks == time) & (adata.obs.Region == region)
        ).astype(str)
        adata.obs.loc[cells,"temp",] = "True"
        adata.obs.loc[adata.obs["temp"] == "True", "temp"] = adata.obs.loc[
            adata.obs["temp"] == "True", "majorclass"
        ]
        adata.obs["temp"] = adata.obs.temp.replace({"False": np.nan})
        sc.pl.umap(adata, color="temp", frameon=False, title="",size = 30)
        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        plt.legend('',frameon=False)
        plt.savefig(
            output_file_path + "HC_"+region + "_" + time + ".svg",
            bbox_inches="tight",
        )
