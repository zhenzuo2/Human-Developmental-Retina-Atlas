# Import packages
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib

matplotlib.rcParams.update({"font.size": 20})

####################################################################################################
cell_type = "AC"
adata_file = "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/AC_annotation_adult/AC_merged_object.h5ad"
meta_cluster_adata_file = "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/AC_annotation_adult/AC_merged_object.h5ad"
subclass_reference_file = "/storage/chentemp/zz4/adult_dev_compare/data/manual_reference/AC_subclass_annotation.csv"
majorclass_reference_file = (
    "/storage/chentemp/zz4/adult_dev_compare/data/cell_subclass_majorclass_mapping.csv"
)

adata = sc.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = sc.read(meta_cluster_adata_file).obs

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.replace(
    dict(zip(subclass_reference.leiden.astype(str), subclass_reference.subclass))
).fillna(cell_type + " Precursor")

adata.obs["meta_cluster_author_cell_type"] = meta_cluster.author_cell_type.astype(
    "category"
)
adata = adata[adata.obs["subclass"] != "Mislabel"]

# Read reference file
majorclass_reference = pd.read_csv(majorclass_reference_file)
adata.obs["majorclass"] = (
    adata.obs["subclass"]
    .map(dict(zip(majorclass_reference.subclass, majorclass_reference.majorclass)))
    .fillna(cell_type + " Precursor")
)

adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")


def find_convertible_ints(input_list):
    convertible_ints = []
    for item in input_list:
        if isinstance(item, str):
            try:
                int_value = int(item)
                convertible_ints.append(item)
            except ValueError:
                pass
    return convertible_ints


adata.obs["subclass"] = adata.obs["subclass"].astype(str)
adata.obs.loc[
    adata.obs["subclass"].isin(find_convertible_ints(set(adata.obs["subclass"]))),
    "subclass",
] = (
    cell_type + " Precursor"
)
matplotlib.rcParams.update({"font.size": 10})
sc.pl.umap(
    adata,
    color="subclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/AC_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
matplotlib.rcParams.update({"font.size": 20})
sc.pl.umap(
    adata,
    color="majorclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/AC_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

mapping_dict = {
    "fetal": "Development",
    "adult": "Adult",
}
adata.obs["sample_source"] = adata.obs["sample_source"].map(mapping_dict)
sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc=None,
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
    palette={"Development": "#1F77B4", "Adult": "#FF7F0E"},
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/AC_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
adata = adata[adata.obs.sample_source == "Development"]
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/filtered_major_class.csv"
)
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index, "Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/AC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
####################################################################################################
cell_type = "BC"
adata_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/"
    + cell_type
    + "_annotation_adult/"
    + cell_type
    + "_merged_object.h5ad"
)
meta_cluster_adata_file = adata_file
subclass_reference_file = (
    "/storage/chentemp/zz4/adult_dev_compare/data/manual_reference/"
    + cell_type
    + "_subclass_annotation.csv"
)

adata = sc.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = sc.read(meta_cluster_adata_file).obs

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.replace(
    dict(zip(subclass_reference.leiden.astype(str), subclass_reference.subclass))
).fillna(cell_type + " Precursor")

adata.obs["meta_cluster_author_cell_type"] = meta_cluster.author_cell_type.astype(
    "category"
)
adata = adata[adata.obs["subclass"] != "Mislabel"]

# Read reference file
majorclass_reference = pd.read_csv(majorclass_reference_file)
adata.obs["majorclass"] = (
    adata.obs["subclass"]
    .map(dict(zip(majorclass_reference.subclass, majorclass_reference.majorclass)))
    .fillna(cell_type + " Precursor")
)

adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")

adata.obs["subclass"] = adata.obs["subclass"].astype(str)
adata.obs.loc[
    adata.obs["subclass"].isin(find_convertible_ints(set(adata.obs["subclass"]))),
    "subclass",
] = (
    cell_type + " Precursor"
)
sc.pl.umap(
    adata,
    color="subclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="majorclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata.obs["sample_source"] = adata.obs["sample_source"].map(mapping_dict)
sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc=None,
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
    palette={"Development": "#1F77B4", "Adult": "#FF7F0E"},
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata = adata[adata.obs.sample_source == "Development"]
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/filtered_major_class.csv"
)
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index, "Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/BC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

####################################################################################################
cell_type = "RGC"
adata_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/"
    + cell_type
    + "_annotation_adult/"
    + cell_type
    + "_merged_object.h5ad"
)
meta_cluster_adata_file = adata_file
subclass_reference_file = (
    "/storage/chentemp/zz4/adult_dev_compare/data/manual_reference/"
    + cell_type
    + "_subclass_annotation.csv"
)

adata = sc.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = sc.read(meta_cluster_adata_file).obs

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.replace(
    dict(zip(subclass_reference.leiden.astype(str), subclass_reference.subclass))
).fillna(cell_type + " Precursor")

adata.obs["meta_cluster_author_cell_type"] = meta_cluster.author_cell_type.astype(
    "category"
)
adata = adata[adata.obs["subclass"] != "Mislabel"]

# Read reference file
majorclass_reference = pd.read_csv(majorclass_reference_file)
adata.obs["majorclass"] = (
    adata.obs["subclass"]
    .map(dict(zip(majorclass_reference.subclass, majorclass_reference.majorclass)))
    .fillna(cell_type + " Precursor")
)

adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")

adata.obs["subclass"] = adata.obs["subclass"].astype(str)
adata.obs.loc[
    adata.obs["subclass"].isin(find_convertible_ints(set(adata.obs["subclass"]))),
    "subclass",
] = (
    cell_type + " Precursor"
)
sc.pl.umap(
    adata,
    color="subclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="majorclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata.obs["sample_source"] = adata.obs["sample_source"].map(mapping_dict)
sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc=None,
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
    palette={"Development": "#1F77B4", "Adult": "#FF7F0E"},
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source == "Development"]
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/filtered_major_class.csv"
)
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index, "Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/RGC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
####################################################################################################
cell_type = "HC"
adata_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/"
    + cell_type
    + "_annotation_adult/"
    + cell_type
    + "_merged_object.h5ad"
)
meta_cluster_adata_file = adata_file
subclass_reference_file = (
    "/storage/chentemp/zz4/adult_dev_compare/data/manual_reference/"
    + cell_type
    + "_subclass_annotation.csv"
)

adata = sc.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = sc.read(meta_cluster_adata_file).obs

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.replace(
    dict(zip(subclass_reference.leiden.astype(str), subclass_reference.subclass))
).fillna(cell_type + " Precursor")

adata.obs["meta_cluster_author_cell_type"] = meta_cluster.author_cell_type.astype(
    "category"
)
adata = adata[adata.obs["subclass"] != "Mislabel"]

# Read reference file
majorclass_reference = pd.read_csv(majorclass_reference_file)
adata.obs["majorclass"] = (
    adata.obs["subclass"]
    .map(dict(zip(majorclass_reference.subclass, majorclass_reference.majorclass)))
    .fillna(cell_type + " Precursor")
)

adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")

adata.obs["subclass"] = adata.obs["subclass"].astype(str)
adata.obs.loc[
    adata.obs["subclass"].isin(find_convertible_ints(set(adata.obs["subclass"]))),
    "subclass",
] = (
    cell_type + " Precursor"
)
sc.pl.umap(
    adata,
    color="subclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="majorclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata.obs["sample_source"] = adata.obs["sample_source"].map(mapping_dict)
sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc=None,
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
    palette={"Development": "#1F77B4", "Adult": "#FF7F0E"},
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source == "Development"]
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/filtered_major_class.csv"
)
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index, "Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/HC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

####################################################################################################
cell_type = "Cone"
adata_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/"
    + cell_type
    + "_annotation_adult/"
    + cell_type
    + "_merged_object.h5ad"
)
meta_cluster_adata_file = adata_file
subclass_reference_file = (
    "/storage/chentemp/zz4/adult_dev_compare/data/manual_reference/"
    + cell_type
    + "_subclass_annotation.csv"
)

adata = sc.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = sc.read(meta_cluster_adata_file).obs

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.replace(
    dict(zip(subclass_reference.leiden.astype(str), subclass_reference.subclass))
).fillna(cell_type + " Precursor")

adata.obs["meta_cluster_author_cell_type"] = meta_cluster.author_cell_type.astype(
    "category"
)
adata = adata[adata.obs["subclass"] != "Mislabel"]

# Read reference file
majorclass_reference = pd.read_csv(majorclass_reference_file)
adata.obs["majorclass"] = (
    adata.obs["subclass"]
    .map(dict(zip(majorclass_reference.subclass, majorclass_reference.majorclass)))
    .fillna(cell_type + " Precursor")
)

adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")

adata.obs["subclass"] = adata.obs["subclass"].astype(str)
adata.obs.loc[
    adata.obs["subclass"].isin(find_convertible_ints(set(adata.obs["subclass"]))),
    "subclass",
] = (
    cell_type + " Precursor"
)
adata = adata[adata.obs.subclass != "Rod"]
sc.pl.umap(
    adata,
    color="subclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="majorclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata.obs["sample_source"] = adata.obs["sample_source"].map(mapping_dict)
sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc=None,
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
    palette={"Development": "#1F77B4", "Adult": "#FF7F0E"},
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source == "Development"]
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/filtered_major_class.csv"
)
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index, "Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/Cone_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


####################################################################################################


cell_type = "Rod"
adata_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/"
    + cell_type
    + "_annotation_adult/"
    + cell_type
    + "_merged_object.h5ad"
)
meta_cluster_adata_file = adata_file
subclass_reference_file = (
    "/storage/chentemp/zz4/adult_dev_compare/data/manual_reference/"
    + cell_type
    + "_subclass_annotation.csv"
)

adata = sc.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = sc.read(meta_cluster_adata_file).obs

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.replace(
    dict(zip(subclass_reference.leiden.astype(str), subclass_reference.subclass))
).fillna(cell_type + " Precursor")

adata.obs["meta_cluster_author_cell_type"] = meta_cluster.author_cell_type.astype(
    "category"
)
adata = adata[adata.obs["subclass"] != "Mislabel"]

# Read reference file
majorclass_reference = pd.read_csv(majorclass_reference_file)
adata.obs["majorclass"] = (
    adata.obs["subclass"]
    .map(dict(zip(majorclass_reference.subclass, majorclass_reference.majorclass)))
    .fillna(cell_type + " Precursor")
)

adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")

adata.obs["subclass"] = adata.obs["subclass"].astype(str)
adata.obs.loc[
    adata.obs["subclass"].isin(find_convertible_ints(set(adata.obs["subclass"]))),
    "subclass",
] = (
    cell_type + " Precursor"
)

adata = adata[adata.obs.subclass != "ML_Cone"]
adata = adata[adata.obs.subclass != "S_Cone"]
sc.pl.umap(
    adata,
    color="subclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="majorclass",
    frameon=False,
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata.obs["sample_source"] = adata.obs["sample_source"].map(mapping_dict)
sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc=None,
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
    palette={"Development": "#1F77B4", "Adult": "#FF7F0E"},
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/"
    + cell_type
    + "_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source == "Development"]
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/filtered_major_class.csv"
)
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index, "Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize=30,
    frameon=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure2/Rod_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
