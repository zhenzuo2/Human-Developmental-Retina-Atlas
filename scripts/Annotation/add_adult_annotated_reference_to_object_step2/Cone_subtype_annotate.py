import scanpy as sc
import matplotlib as mpl
import pandas as pd

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/Cone_subtype.csv"
)
df["celltype"] = df["subclass"]
df["subclass"] = df["majorclass"]
df["majorclass"] = "Cone"

Cone = adata[df["Unnamed: 0"]]

Cone.obs["celltype"] = df["subclass"].values
Cone.obs["subclass"] = df["subclass"].values
Cone.obs["celltype"] = df["celltype"].values

Cone.obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/Cone_subtype_processed.csv"
)
