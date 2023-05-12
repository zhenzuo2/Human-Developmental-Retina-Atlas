import sys
import scvelo as scv
import scanpy as sc
import scvi
import anndata
import pandas as pd
import os
import tempfile
import anndata as ad

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_AC.csv"
)
nrpc = adata[meta.loc[meta.majorclass=="NRPC","Unnamed: 0"]]
nrpc.obs['subclass'] = "NRPC"
subtypes = scv.read("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/AC_subtype.h5ad")

adata = ad.concat([subtypes, nrpc])

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
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp

adata = run_umap_scvi(adata)

adata.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/AC_subtype_NRPC.h5ad"
)
adata.obs.to_csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/umap/AC_subtype_NRPC_obs.csv")
