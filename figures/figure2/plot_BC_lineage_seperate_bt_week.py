# Import packages
# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns

sc.set_figure_params(transparent=True, fontsize=15)
NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC = NRPC[NRPC.obs.leiden!="12"]
NRPC.obs["subclass"] = "Unknown"

sc.pp.neighbors(NRPC, use_rep="X_scVI")
sc.tl.leiden(NRPC, resolution=3)

clusters = ["28", "3", "13", "15"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["6", "30","35","25"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["21","14","27"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["33","0"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["12","22"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["4", "19","23"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"


adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_BC/adata_umap.h5ad")
adata.obs["subclass"] = adata.obs["majorclass"].astype(str)
cells = [x for x in NRPC[NRPC.obs.leiden == "14"].obs.index if x in adata.obs.index]
adata.obs.loc[cells,'subclass'] = "Cluster 2"
cells = [x for x in NRPC[NRPC.obs.leiden == "27"].obs.index if x in adata.obs.index]
adata.obs.loc[cells,'subclass'] = "Cluster 2"
cells = [x for x in NRPC[NRPC.obs.leiden == "21"].obs.index if x in adata.obs.index]
adata.obs.loc[cells,'subclass'] = "Cluster 1"


df = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/BC_subtype.csv")
df.index = df["Unnamed: 0"].values
df = df.loc[[x for x in df.index if x in adata.obs.index],:]
adata.obs.loc[df.index, "subclass"] = df.majorclass


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
#cells = []
#for sp in set(adata.obs.subclass):
#    cells = cells + [adata[adata.obs.subclass == sp].obs.index[0]]

adata.obs["subclass"] = adata.obs["subclass"].replace({'ON-BC':'ON Cone BC', 'OFF-BC':'OFF Cone BC'})
for region in set(adata.obs.Region):
    for weeks in set(adata.obs.Weeks):
        adata.obs["temp"] = np.nan
        subset = (adata.obs.Region == region) & (adata.obs.Weeks == weeks)
        #adata.obs.loc[cells, "temp"] = adata.obs.loc[cells, "subclass"]
        adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "subclass"]
        adata.obs["temp"] = pd.Categorical(
            list(adata.obs["temp"]),
            categories=[
                'BC Precursor', 'OFF Cone BC', 'ON Cone BC', 'RBC',"Cluster 1","Cluster 2"
            ],
        )
        sc.pl.umap(
            adata,
            color="temp",
            size=80,
            title="",
            frameon=False,
            legend_loc="None",
            palette={
                "BC Precursor": "#1f77b4",
                "OFF Cone BC": "#2ca02c",
                "ON Cone BC": "#d62728",
                "RBC": "#9467bd",
                "Cluster 2": "#8c564b",
                "Cluster 1": "#e377c2",
            },
        )
        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        plt.savefig(
            "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/BC_lineage_umap_by_cell_type_"
            + region
            + "_"
            + weeks
            + ".png",
            dpi=600,
            bbox_inches="tight",
            transparent=True,
            #backend="cairo",
        )
