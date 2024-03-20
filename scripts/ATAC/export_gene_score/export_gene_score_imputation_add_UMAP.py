import numpy as np
import pandas as pd
import scanpy as sc

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_score_imputed.h5ad"
)
df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/ATAC_umap/ATAC_UMAPHarmony_cor.csv"
)
df.index = df["Unnamed: 0"].values
df.index = [x.replace("#", "_") for x in df.index.values]

df = df.loc[adata.obs.index,:]
x = df["Harmony#UMAP_Dimension_1"]
y = df["Harmony#UMAP_Dimension_2"]

adata.obsm["X_umap"] = np.array(list(zip(x, y)))
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_score_imputed_umap.h5ad"
)
