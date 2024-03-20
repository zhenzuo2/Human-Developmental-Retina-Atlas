import scanpy as sc
import matplotlib as mpl
import pandas as pd

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/HC_subtype.csv"
)
df["celltype"] = df["subclass"]
df["subclass"] = df["majorclass"]
df["majorclass"] = "HC"

HC = adata[df["Unnamed: 0"]]

HC.obs["celltype"] = df["subclass"].values
HC.obs["subclass"] = df["subclass"].values
HC.obs["celltype"] = df["celltype"].values

HC.obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/HC_subtype_processed.csv"
)
