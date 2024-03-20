import scanpy as sc
import pandas as pd
import tempfile
import scvi
import scvelo as scv
import os
import sys
x = sys.argv[1]

df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/data/cell_cycle_genes.csv")
cell_cycle_genes = list(df.Gene)

# Read data
def run_umap_scvi(adata, batch_key="sampleid", labels_key="subclass"):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata = adata.copy()
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
    try:
        scvi.settings.seed = 0
        scvi.model.SCVI.setup_anndata(adata, labels_key=labels_key, batch_key=batch_key)
        vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
        vae.train()
        adata.obsm["X_scVI"] = vae.get_latent_representation()
    except:
        scvi.settings.seed = 0
        scvi.model.SCVI.setup_anndata(adata, labels_key=labels_key)
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

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated.h5ad"
)
if x != "ALL":
    cell1 = list(adata.obs.index.values)
    df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/"+x+".csv")
    cell2 = list(df['Unnamed: 0'])
    cells = list(set(cell1) & set(cell2))
    adata = adata[cells,:]
    if x =="RPC_MG":
        shared_genes = list(set(adata.var.index) - set(cell_cycle_genes))
        adata = adata[:, shared_genes]
adata
adata = run_umap_scvi(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/"+x+".h5ad"
)
