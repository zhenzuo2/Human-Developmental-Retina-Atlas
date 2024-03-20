# Import packages
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import matplotlib
import seaborn as sns

matplotlib.font_manager._load_fontmanager(try_read_cache=False)
plt.rcParams["font.family"] = "Arial"

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC_w_NRPC.h5ad"
)
adata = adata[adata.obs.majorclass == "NRPC"]
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

adata.obs["Weeks"] = adata.obs.Days.map(
    {
        59: "PCW8",
        70: "PCW10",
        76: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW15",
        103: "PCW15",
        116: "PCW15",
        137: "PCW19",
        141: "PCW19",
        142: "PCW19",
        162: "PCW23",
        165: "PCW23",
    }
)
adata.obs["ONECUT1"] = adata[:, "ONECUT1"].X.toarray()
adata = adata[
    adata.obs.Weeks.isin(
        [
            "PCW8",
            "PCW10",
            "PCW13",
            "PCW15",
            "PCW19",
        ]
    )
]
sns.violinplot(
    data=adata.obs,
    x="Weeks",
    y="ONECUT1",
    hue="Weeks",
    inner=None,
    dodge = False,
    order = [
            "PCW8",
            "PCW10",
            "PCW13",
            "PCW15",
            "PCW19",
        ],
    palette={
        "PCW8": "#FF0000",
        "PCW10": "#FF4500",
        "PCW13": "#FFA500",
        "PCW15": "#4169E1",
        "PCW19": "#0000FF",
    },
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
# plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
plt.legend("", frameon=False)
plt.xticks(rotation="vertical")
plt.xlabel("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/sns.violinplot_ONECUT1.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
plt.clf()
