import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scvelo as scv

adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/AC.h5ad")
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/"

meta_file = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/AC_subtype.csv")
meta_file.index = meta_file["Unnamed: 0"]
common_cells = [x for x in meta_file.index if x in adata.obs.index]
meta_file = meta_file.loc[common_cells, :]
adata.obs["subclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
for sublcass in set(meta_file.subclass):
    adata.obs.loc[
        meta_file.loc[meta_file.subclass == sublcass, :].index, "subclass"
    ] = sublcass
    
for majorclass in set(meta_file.majorclass):
    adata.obs.loc[
        meta_file.loc[meta_file.majorclass == majorclass, :].index, "majorclass"
    ] = majorclass

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
for sp in set(adata.obs.majorclass):
    cells = cells + [adata[adata.obs.majorclass == sp].obs.index[0]]

for time in set(adata.obs["Weeks"]):
    for region in set(adata.obs.Region):
        adata.obs["temp"] = (
            (adata.obs.Weeks == time) & (adata.obs.Region == region)
        ).astype(str)
        adata.obs.loc[
            cells,
            "temp",
        ] = "True"
        adata.obs.loc[adata.obs["temp"] == "True", "temp"] = adata.obs.loc[
            adata.obs["temp"] == "True", "majorclass"
        ]
        adata.obs["temp"] = adata.obs.temp.replace({"False": np.nan})
        sc.pl.umap(adata, color="temp", frameon=False, title="",size = 60.0)
        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        plt.legend('',frameon=False)
        plt.savefig(
            output_file_path + region + "_" + time + ".svg",
            bbox_inches="tight",
        )
