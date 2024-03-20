import scanpy as sc
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import pandas as pd

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

df = pd.DataFrame(adata.obs)
df.loc[:, "cell_id"] = list(df.index.values)
grouped_data = df.groupby("majorclass")
# Define the number of rows to downsample
num_rows = 4000
# Sample the same number of rows from each group
downsampled_data = grouped_data.apply(
    lambda x: x.sample(n=num_rows, random_state=42, replace=True)
)
# Reset the index
downsampled_data = downsampled_data.reset_index(drop=True)

adata_subset = adata[downsampled_data.cell_id]
adata_subset.obs["majorclass"] = pd.Categorical(
    list(adata_subset.obs["majorclass"]),
    categories=["BC", "Cone", "Rod", "MG", "NRPC", "PRPC", "RGC", "AC", "HC"],
)

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

df = pd.DataFrame(adata.obs)
df.loc[:, "cell_id"] = list(df.index.values)
grouped_data = df.groupby("majorclass")
# Define the number of rows to downsample
num_rows = 4000
# Sample the same number of rows from each group
downsampled_data = grouped_data.apply(
    lambda x: x.sample(n=num_rows, random_state=42, replace=True)
)
# Reset the index
downsampled_data = downsampled_data.reset_index(drop=True)

adata_subset = adata[downsampled_data.cell_id]

marker_genes_dict = {
    "OCA": ["TYR", "OCA2", "TYRP1", "SLC45A2", "SLC24A5", "LRMDA"],
    "OA": "GPR143",
    "HPS": [
        "HPS1",
        "AP3B1",
        "HPS3",
        "HPS4",
        "HPS5",
        "HPS6",
        "DTNBP1",
        "BLOC1S3",
        "BLOC1S6",
        "AP3D1",
        "BLOC1S5",
    ],
    "CHS": "LYST",
    "FHONDA": "SLC38A8",
    "Aniridia": "PAX6",
    "FRMD7-related infantile nystagmus": "FRMD7",
    "AHR-related FH and infantile nystagmus": "AHR",
    "Achromatopsia": ["CNGB3", "CNGA3", "GNAT2", "PDE6C", "PDE6H", "ATF6"],
}

sc.pl.heatmap(
    adata_subset,
    marker_genes_dict,
    groupby="majorclass",
    var_group_labels=True,
    dendrogram=True,
    vmax=1.5,
)
plt.xlabel("")
plt.ylabel("")
plt.title("")
fig = plt.gcf()
fig.set_size_inches(15, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure7/FovealHypoplasiaHetamap.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
