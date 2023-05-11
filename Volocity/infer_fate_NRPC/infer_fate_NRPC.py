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
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, batch_key="sampleid")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp


adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered.h5ad"
)
meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv"
)

meta = meta[(meta.majorclass == "NRPC") & (meta.Time != "Adult")]
adata = adata[meta["Unnamed: 0"]]

adata = run_umap_scvi(adata)

fate = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/to_terminal_states.csv"
)
fate.index = fate["Unnamed: 0"].values

fate = fate.drop(columns = 'Unnamed: 0')
cell_types = fate.columns.values

adata = adata[[x for x in fate.index if x in adata.obs.index]]
fate = fate.loc[list(adata.obs.index.values), cell_types]

for x in cell_types:
    adata.obs[x] = fate[x]

adata.obs["terminal_state"] = adata.obs.loc[:, cell_types].idxmax(axis=1)
adata.obs["terminal_state_p"] = adata.obs.loc[:, cell_types].max(axis=1)

adata.obs.to_csv("/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.csv")
adata.write("/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad")