import pandas as pd
import scvelo as scv
import scanpy as sc

AC = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/AC.csv"
)
BC = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/BC.csv"
)
Rod = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/Rod.csv"
)
Cone = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/Cone.csv"
)
HC = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/HC.csv"
)
RGC = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/RGC.csv"
)
MG = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/MG.csv"
)

meta = pd.concat([AC, BC, Rod, Cone, HC, RGC, MG])
meta.index = meta["Unnamed: 0"].values
meta = meta.drop(columns=["Unnamed: 0", "meta_cluster_author_cell_type", "leiden"])

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000.h5ad"
)

adata = adata[meta.index.values]
adata.obs = meta

adata.write(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad"
)
meta.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.obs.csv"
)

meta.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/merged_raw_filtered_umap_10000_major_sub_class.obs.csv"
)