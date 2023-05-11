import scvelo as scv
import pandas as pd

cell_type = "BC"
adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/BC.h5ad"
meta_cluster_adata_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/BC_annotation_adult/BC_merged_object.h5ad"
subclass_reference_file = "/storage/singlecell/zz4/fetal_snakemake/data/manual_reference/BC_subclass_annotation.csv"
majorclass_reference_file = (
    "/storage/singlecell/zz4/fetal_snakemake/data/cell_subclass_majorclass_mapping.csv"
)

adata = scv.read(adata_file)
# Read meta_cluster file with cluster information
meta_cluster = scv.read(meta_cluster_adata_file).obs.loc[adata.obs.index, :]
meta_cluster["author_cell_type"] = meta_cluster.author_cell_type.astype(str).astype(int)

# Read reference file
subclass_reference = pd.read_csv(subclass_reference_file)

adata.obs["subclass"] = meta_cluster.author_cell_type.map(
    dict(zip(subclass_reference.leiden, subclass_reference.subclass))
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
adata.obs.to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/BC_subtype.csv"
)
adata.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/BC_subtype.h5ad"
)
