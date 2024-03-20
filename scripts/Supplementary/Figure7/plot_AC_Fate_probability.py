import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC_w_NRPC.h5ad")
df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/AC_infer_fate/AC_to_terminal_states.csv")
df.index = df["Unnamed: 0"].values

adata.obs["SACs"] = df.loc[adata.obs.index,"SACs"]
adata.obs["GABAergic"] = df.loc[adata.obs.index,"GABAergic"]
adata.obs["Glycinergic"] = df.loc[adata.obs.index,"Glycinergic"]

sc.pl.umap(adata[adata.obs.subclass.isin(["NRPC",])],color ="Glycinergic",vmax = 0.00001,frameon = False, title = "Glycinergic Fate Probabilities")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure7/Glycinergic_fate_p.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)

sc.pl.umap(adata[adata.obs.subclass.isin(['AC Precursor',])],color ="Glycinergic",vmax = 1,frameon = False, title = "Glycinergic Fate Probabilities")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure7/Glycinergic_fate_p_Precursor.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)

sc.pl.umap(adata[adata.obs.subclass.isin(["NRPC",])],color ="GABAergic",vmax = 0.016,frameon = False,title ="GABAergic Fate Probabilities")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure7/GABAergic_fate_p.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)

sc.pl.umap(adata[adata.obs.subclass.isin(["AC Precursor",])],color ="GABAergic",vmax = 1,frameon = False,title ="GABAergic Fate Probabilities")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure7/GABAergic_fate_p_Precursor.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)