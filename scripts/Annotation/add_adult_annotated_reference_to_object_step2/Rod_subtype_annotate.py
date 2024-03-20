import scanpy as sc
import matplotlib as mpl
import pandas as pd

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/Rod_subtype.csv"
)
df["celltype"] = df["subclass"]
df["subclass"] = df["majorclass"]
df["majorclass"] = "Rod"

Rod = adata[df["Unnamed: 0"]]

Rod.obs["celltype"] = df["subclass"].values
Rod.obs["subclass"] = df["subclass"].values
Rod.obs["celltype"] = df["celltype"].values

Rod.obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/Rod_subtype_processed.csv"
)
