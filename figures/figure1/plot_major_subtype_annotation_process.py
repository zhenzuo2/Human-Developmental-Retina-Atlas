# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_wadult_umap_10000.h5ad"
)
scv.pl.umap(
    adata, color="batch", size=50, legend_loc="right margin", title="", sort_order=False
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_with_adult_by_batch.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

scv.pl.umap(
    adata,
    color="majorclass",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_with_adult_by_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
####################################################################################################
cell_type = "AC"
adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/AC_annotation_adult/AC_merged_object.h5ad"
meta_cluster_adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/AC_annotation_adult/AC_merged_object.h5ad"
subclass_reference_file = "/storage/singlecell/zz4/fetal_snakemake/data/manual_reference/AC_subclass_annotation.csv"
majorclass_reference_file = (
    "/storage/singlecell/zz4/fetal_snakemake/data/cell_subclass_majorclass_mapping.csv"
)

adata = scv.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs

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
scv.pl.umap(
    adata,
    color="subclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20"
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/AC_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

scv.pl.umap(
    adata,
    color="majorclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/AC_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
    palette = {"fetal": "#7F7F7F", "adult": "#EF553B"}
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/AC_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
adata = adata[adata.obs.sample_source=="fetal"]
df = pd.read_csv('/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv')
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index,"Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/AC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
####################################################################################################
cell_type = "BC"
adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/"+cell_type+"_annotation_adult/"+cell_type+"_merged_object.h5ad"
meta_cluster_adata_file = adata_file
subclass_reference_file = "/storage/singlecell/zz4/fetal_snakemake/data/manual_reference/"+cell_type+"_subclass_annotation.csv"

adata = scv.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs

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
scv.pl.umap(
    adata,
    color="subclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20"
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

scv.pl.umap(
    adata,
    color="majorclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
    palette = {"fetal": "#7F7F7F", "adult": "#EF553B"}
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata = adata[adata.obs.sample_source=="fetal"]
df = pd.read_csv('/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv')
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index,"Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/BC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

####################################################################################################
cell_type = "RGC"
adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/"+cell_type+"_annotation_adult/"+cell_type+"_merged_object.h5ad"
meta_cluster_adata_file = adata_file
subclass_reference_file = "/storage/singlecell/zz4/fetal_snakemake/data/manual_reference/"+cell_type+"_subclass_annotation.csv"

adata = scv.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs

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
scv.pl.umap(
    adata,
    color="subclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20"
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

scv.pl.umap(
    adata,
    color="majorclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
    palette = {"fetal": "#7F7F7F", "adult": "#EF553B"}
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source=="fetal"]
df = pd.read_csv('/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv')
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index,"Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/RGC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
####################################################################################################
cell_type = "HC"
adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/"+cell_type+"_annotation_adult/"+cell_type+"_merged_object.h5ad"
meta_cluster_adata_file = adata_file
subclass_reference_file = "/storage/singlecell/zz4/fetal_snakemake/data/manual_reference/"+cell_type+"_subclass_annotation.csv"

adata = scv.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs

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
scv.pl.umap(
    adata,
    color="subclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20"
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

scv.pl.umap(
    adata,
    color="majorclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
    palette = {"fetal": "#7F7F7F", "adult": "#EF553B"}
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source=="fetal"]
df = pd.read_csv('/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv')
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index,"Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/HC_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

####################################################################################################
cell_type = "Cone"
adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/"+cell_type+"_annotation_adult/"+cell_type+"_merged_object.h5ad"
meta_cluster_adata_file = adata_file
subclass_reference_file = "/storage/singlecell/zz4/fetal_snakemake/data/manual_reference/"+cell_type+"_subclass_annotation.csv"

adata = scv.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs

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
scv.pl.umap(
    adata,
    color="subclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20"
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

scv.pl.umap(
    adata,
    color="majorclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
    palette = {"fetal": "#7F7F7F", "adult": "#EF553B"}
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source=="fetal"]
df = pd.read_csv('/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv')
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index,"Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/Cone_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


####################################################################################################


cell_type = "Rod"
adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/"+cell_type+"_annotation_adult/"+cell_type+"_merged_object.h5ad"
meta_cluster_adata_file = adata_file
subclass_reference_file = "/storage/singlecell/zz4/fetal_snakemake/data/manual_reference/"+cell_type+"_subclass_annotation.csv"

adata = scv.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs

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

adata = adata[adata.obs.subclass != 'ML_Cone']
adata = adata[adata.obs.subclass != 'S_Cone']
scv.pl.umap(
    adata,
    color="subclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20"
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_subclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

scv.pl.umap(
    adata,
    color="majorclass",
    size=50,
    legend_loc="on data",
    title="",
    sort_order=False,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_author_majorclass.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="sample_source",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
    palette = {"fetal": "#7F7F7F", "adult": "#EF553B"}
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/"+cell_type+"_sample_source.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


adata = adata[adata.obs.sample_source=="fetal"]
df = pd.read_csv('/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv')
df.index = df["Unnamed: 0"].values
adata.obs["Days"] = df.loc[adata.obs.index,"Days"].values
sc.pl.umap(
    adata,
    color="Days",
    size=50,
    legend_loc="right margin",
    title="",
    sort_order=False,
    legend_fontsize = 30,
    frameon = False,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/Rod_Days.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
