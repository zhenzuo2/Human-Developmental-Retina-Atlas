import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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
for time in set(adata.obs["Weeks"]):
    for region in set(adata.obs.Region):
        adata.obs["temp"] = (
            (adata.obs.Weeks == time) & (adata.obs.Region == region)
        ).astype(str)
        adata.obs.loc[
            [
                "Multi_Fetal_11w2d_FR_2_ACGATTCAGGCTATGT-1",
                "Multi_Fetal_11w2d_NR_AAAGCTTGTCTAACCT-1",
                "Multi_Fetal_11w2d_FR_2_GACGCAACAATATAGG-1",
                "Multi_Fetal_13W_NR_AGCTAAACAATATAGG-1",
                "Multi_Fetal_19W4d_FR_AAACAGCCAATAACCT-1",
                "Multi_Fetal_19W4d_NR_AGTTGCAGTGGGTGAA-1",
                "Multi_Fetal_11w2d_FR_2_TTTAACCTCAAATCGC-1",
                "Multi_Fetal_11w2d_NR_AAACCGAAGGAGTCGG-1",
                "Multi_Fetal_13W_FR_AAATGGCCAGTTTGGC-1",
                "Multi_Fetal_13W_NR_AACAGGATCGATTTGA-1",
                "Multi_Fetal_19W4d_FR_AAACCAACACCAGGTT-1",
                "Multi_Fetal_19W4d_NR_AGAACCAAGGGTTAGA-1",
            ],
            "temp",
        ] = "True"
        adata.obs.loc[adata.obs["temp"] == "True", "temp"] = adata.obs.loc[
            adata.obs["temp"] == "True", "majorclass"
        ]
        adata.obs["temp"] = adata.obs.temp.replace({"False": np.nan})
        sc.pl.umap(adata, color="temp", frameon=False, title="",size = 30)
        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        plt.legend('',frameon=False)
        plt.savefig(
            output_file_path + region + "_" + time + ".svg",
            bbox_inches="tight",
        )
