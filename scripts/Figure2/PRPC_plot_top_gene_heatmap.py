import joblib
import scvelo as scv
import plotly.express as px
import scanpy as sc
import hotspot
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib
from scipy.cluster.hierarchy import leaves_list
matplotlib.rcParams.update({'font.size': 12})

output_file_path = "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/"
hs = joblib.load(
    "/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/PRPC_hs.create_modules_velocity_pseudotime.pkl"
)
adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

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

ii = leaves_list(hs.linkage)

mod_reordered = hs.modules.iloc[ii]

mod_map = {}
y = np.arange(modules.size)

for x in mod_reordered.unique():
    if x == -1:
        continue

mod_map[x] = y[mod_reordered == x].mean()

adata_result = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/PRPC_hs.adata_velocity_pseudotime.h5ad")

sc.pp.normalize_total(adata_result)
sc.pp.scale(adata_result)

top_genes = mod_reordered[mod_reordered==3][:20].index
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='velocity_pseudotime',color_map = "RdBu_r")
fig = plt.gcf()
fig.set_size_inches(10,4)
plt.savefig(
    output_file_path + "PRPC_gene_module_3_heatmap.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi = 300
)

top_genes = mod_reordered[mod_reordered==2][:20].index
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='velocity_pseudotime',color_map = "RdBu_r")
fig = plt.gcf()
fig.set_size_inches(10,4)
plt.yticks(rotation = 0)
plt.savefig(
    output_file_path + "PRPC_gene_module_2_heatmap.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi = 300
)

top_genes = ["C21orf58", "MIS18BP1", "LINC01572", "SMC4", "RTKN2", "ASPM", "KIF14", "APOLD1", "MKI67", "TOP2A"] + ["LINC01414", "PCDH7", "NALCN-AS1", "KIAA1217", "ROBO2", "CNTN5", "MEG8", "FIGN", "PLEKHG1", "ESRRG"]
scv.pl.heatmap(adata_result, var_names=top_genes, sortby='velocity_pseudotime',color_map = "RdBu_r")
fig = plt.gcf()
fig.set_size_inches(10,4)
plt.savefig(
    output_file_path + "PRPC_gene_module_1_heatmap.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi = 300
)