import matplotlib.pyplot as plt
import scvelo as scv
import scanpy as sc
import seaborn as sns
import pandas as pd

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap/merged_h5ad_adult_annotated_umap_X_scANVI.h5ad",
    cache=False,
)

adata = adata[adata.obs.majorclass == "NRPC"]

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=True)
scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)
sc.tl.umap(adata)

fate = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/to_terminal_states.csv"
)
fate.index = fate["Unnamed: 0"].values
adata = adata[[x for x in adata.obs.index if x in fate.index]]
fate = fate.loc[
    list(adata.obs.index.values), ["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]
]

for x in ["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]:
    adata.obs[x] = fate[x]

adata.obs["terminal_state"] = adata.obs.loc[
    :, ["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]
].idxmax(axis=1)
adata.obs["terminal_state_p"] = adata.obs.loc[
    :, ["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]
].max(axis=1)

adata.obs.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/NRPC_fate/NRPC_fate.csv"
)

adata.write(
    "/storage/singlecell/zz4/fetal_bash/results/NRPC_fate/NRPC_fate.h5ad"
)