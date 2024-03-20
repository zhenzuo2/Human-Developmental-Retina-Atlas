import scvelo as scv
import pandas as pd

cell_type = "Rod"
adata_file = "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
meta_cluster_adata_file = "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/Rod_annotation_adult/Rod_merged_object.h5ad"
subclass_reference_file = "/storage/chentemp/zz4/adult_dev_compare/data/manual_reference/Rod_subclass_annotation.csv"
majorclass_reference_file = (
    "/storage/chentemp/zz4/adult_dev_compare/data/cell_subclass_majorclass_mapping.csv"
)

adata = scv.read(adata_file)
adata = adata[adata.obs.majorclass==cell_type]
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
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/Rod_subtype.csv"
)

