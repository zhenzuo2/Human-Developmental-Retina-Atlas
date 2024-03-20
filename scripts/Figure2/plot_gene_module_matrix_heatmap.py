import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib
from scipy.cluster.hierarchy import leaves_list
import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 30
adata_result = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

hs = joblib.load(
    "/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/PRPC_hs.create_modules_velocity_pseudotime.pkl"
)
adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/PRPC_hs.adata_velocity_pseudotime.h5ad"
)

adata.obsm["X_umap"] = adata_result[adata.obs.index].obsm["X_umap"]

modules = hs.create_modules(min_gene_threshold=160, core_only=True, fdr_threshold=0.01)

hs.modules = hs.modules.replace({3: 1, 4: 2, 2: 3})
for x in [
    "GNB4",
    "TEAD4",
    "HAUS2",
    "VKORC1L1",
    "LMNB1-DT",
    "PCGF6",
    "EIF2AK3",
    "FAM185A",
    "AC122719.3",
    "NEDD4",
    "LRRC8B",
    "SPDYA",
    "AC092910.3",
    "SV2C",
    "SOX11",
    "ABHD17B",
]:
    hs.modules[x] = 1
hs.local_correlation_z.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_plot_local_correlations.csv")
hs.plot_local_correlations()
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_plot_local_correlations.tiff",
    dpi=600,
    transparent=True,
    bbox_inches="tight",
)

ii = leaves_list(hs.linkage)

mod_reordered = hs.modules.iloc[ii]

mod_map = {}
y = np.arange(modules.size)

for x in mod_reordered.unique():
    if x == -1:
        continue

mod_map[x] = y[mod_reordered == x].mean()

mod_reordered.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/mod_reordered.csv"
)

module_scores = hs.calculate_module_scores()

adata.obs["Module1"] = module_scores.loc[adata.obs.index, 1].values
adata.obs["Module2"] = module_scores.loc[adata.obs.index, 2].values
adata.obs["Module3"] = module_scores.loc[adata.obs.index, 3].values

adata.obs.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_Module_gene_score.csv")

plt.rcParams["font.size"] = 10

plt.clf()
sc.pl.umap(adata, color="Module1", vmax=6, frameon=False)
fig = plt.gcf()
plt.ylabel("")
plt.xlabel("")
plt.title("Module 1", fontsize=40)
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_Module1_gene_score.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)

plt.clf()
sc.pl.umap(adata, color="Module2", vmax=2, frameon=False)
fig = plt.gcf()
plt.ylabel("")
plt.xlabel("")
plt.title("Module 2", fontsize=40)
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_Module2_gene_score.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)

plt.clf()
sc.pl.umap(adata, color="Module3", vmax=2, frameon=False)
fig = plt.gcf()
plt.ylabel("")
plt.xlabel("")
plt.title("Module 3", fontsize=40)
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_Module3_gene_score.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)
