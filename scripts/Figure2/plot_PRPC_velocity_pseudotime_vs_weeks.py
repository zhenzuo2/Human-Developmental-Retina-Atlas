import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

font = {"family": "Arial", "size": 24}
plt.rc("font", **font)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/infer_velocity_pseudotime/PRPC.csv"
)

df["Days"] = df["Days"].astype(float)
df["Weeks"] = df.Days.map(
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
df["Region"] = df["Region"].astype("category")
df["Region"].cat.reorder_categories(["Whole Eye", "Macula", "Peripheral"])
df.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/plot_PRPC_velocity_pseudotime_vs_weeks.csv")
my_pal = {"Whole Eye": "#ffdd05", "Macula" : "#F8766D", "Peripheral": "#00BFC4"}
sns.boxplot(
    x="Weeks",
    y="velocity_pseudotime",
    hue="Region",
    data=df,
    order=["PCW8", "PCW10", "PCW13", "PCW15", "PCW19", "PCW23"],
    palette=my_pal,
    showfliers=False,
)
plt.ylabel("Gene-shared latent time")
plt.xlabel("Post Conception Week")
plt.xticks(rotation=45)
plt.legend()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/plot_PRPC_velocity_pseudotime_vs_weeks.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata = adata[df["Unnamed: 0"]]
adata.obs["velocity_pseudotime"] = list(df["velocity_pseudotime"].astype(float))
font = {"family": "Arial", "size": 20}
plt.rc("font", **font)
sc.pl.umap(
    adata,
    color="velocity_pseudotime",
    color_map="plasma",
    s=10,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.title("")
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/plot_PRPC_velocity_pseudotime_UMAP.tiff",
    transparent=True,
    dpi=300,
)
