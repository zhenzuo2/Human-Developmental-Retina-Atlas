import numpy as np
import scanpy as sc
import pandas as pd
import magic
import re
import matplotlib.pyplot as plt

def is_valid_string(s):
    return not s.startswith("MT") and re.match("^[a-zA-Z0-9_]*$", s) is not None

def smooth(adata):
    magic_op = magic.MAGIC()
    magic_op.set_params(n_jobs=10)
    emt_magic = magic_op.fit_transform(adata.X, genes="all_genes")
    emt_magic = magic_op.transform(genes="all_genes")
    adata.X = emt_magic
    return adata

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad")
adata.obs["Region"] = adata.obs["Region"].replace({"Peripheral":"Periphery"})
adata = adata[(adata.obs.majorclass == "PRPC") & (adata.obs.Region != "Whole Eye")]
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=False)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
hvg = list(adata.var.loc[adata.var.highly_variable].index)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/PRPC_monocle3_DE_analysis/compare_mod_region.csv"
)
df = df.loc[(df.q_value < 0.01) & (df.num_cells_expressed > 5000), :]
gene_list_1 = list(df.gene_short_name)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/PRPC_monocle3_DE_analysis/fit_coefs_region_models.csv"
)
df = df.loc[
    (df.term == "RegionPeripheral")
    & (df.q_value < 0.01)
    & (df.normalized_effect > 2.15),
]
gene_list_2 = list(df.gene_short_name)
gene_list = list(set(gene_list_1).intersection(gene_list_2))

filtered_strings = [s for s in gene_list if is_valid_string(s)]
filtered_strings = [s for s in filtered_strings if s in hvg]
filtered_strings = [s for s in filtered_strings if s != "WIF1"]
# Print the result
print(len(filtered_strings))
print(filtered_strings)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/PRPC_monocle3_DE_analysis/compare_mod_region.csv"
)
df = df.loc[(df.q_value < 0.01) & (df.num_cells_expressed > 5000), :]
gene_list_1 = list(df.gene_short_name)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/PRPC_monocle3_DE_analysis/fit_coefs_region_models.csv"
)
df = df.loc[
    (df.term == "RegionPeripheral")
    & (df.q_value < 0.01)
    & (df.normalized_effect < -1.08),
]
gene_list_2 = list(df.gene_short_name)
gene_list = list(set(gene_list_1).intersection(gene_list_2))

# Filter the list
filtered_strings2 = [s for s in gene_list if is_valid_string(s)]
filtered_strings2 = [s for s in filtered_strings2 if s in hvg]
# Print the result
print(len(filtered_strings2))
print(filtered_strings2)

markers = {"Macula": filtered_strings2, "Periphery": filtered_strings}
all_values = markers["Macula"] + markers["Periphery"]

adata = adata[:, all_values]
adata = smooth(adata)

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

adata.obs["Region_Weeks"] = pd.concat(
    [adata.obs["Region"].astype(str), adata.obs["Weeks"]], axis=1
).agg("_".join, axis=1)

samples_per_group = 2000
# Downsample the DataFrame
downsampled_df = adata.obs.groupby("Region_Weeks", group_keys=False).apply(
    lambda group: group.sample(samples_per_group, replace=True)
)

sc.pl.heatmap(
    adata[downsampled_df.index],
    markers,
    groupby="Region_Weeks",
    dendrogram=False,
    swap_axes=False,
    vmax=3,
)
plt.ylabel("")
plt.xlabel("")
plt.title("")
plt.xticks(fontsize=30)
fig = plt.gcf()
fig.set_size_inches(10, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/plot_PRPC_Deg_heatmap.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
