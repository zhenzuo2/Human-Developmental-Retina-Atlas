import matplotlib.pyplot as plt
import scvelo as scv
import scanpy as sc
import seaborn as sns
import pandas as pd
import tempfile
import scvi
import os


def run_umap_scvi(adata):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=2000, subset=True
    )
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=2)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp


adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap/merged_h5ad_adult_annotated_umap_X_scANVI.h5ad",
    cache=False,
)

adata = adata[adata.obs.majorclass == "NRPC"]

adata = run_umap_scvi(adata)

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

adata.obs.to_csv("/storage/singlecell/zz4/fetal_bash/results/NRPC_fate/NRPC_fate.csv")

adata.write("/storage/singlecell/zz4/fetal_bash/results/NRPC_fate/NRPC_fate.h5ad")