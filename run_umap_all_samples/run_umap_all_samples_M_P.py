import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scvi
import scanpy as sc
import tempfile
import scvelo as scv
import os


sc.set_figure_params(figsize=(10, 10))
scvi.settings.seed = 0

adata = sc.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw.h5ad"
)

meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/merged_raw_filtered_meta_major_class.csv"
)
meta.index = meta["Unnamed: 0"].values

adata = adata[
    meta["Unnamed: 0"].values,
]

adata.obs["majorclass"] = meta.loc[adata.obs.index.values, "majorclass"]

adata = adata[
    adata.obs.majorclass.isin(["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]),
]

def run_umap_scvi(adata):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, batch_key="sampleid", labels_key="majorclass")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    scvi.settings.seed = 0

    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key="majorclass",
        unlabeled_category="Unknown",
    )
    lvae.train(max_epochs=100, batch_size=640)
    adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp

adata = adata[adata.obs.Region.isin(['Macula','Peripheral'])&(adata.obs['Data Type']=='Multiomics')]
adata = run_umap_scvi(adata)
adata.write(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_M_P.h5ad"
)