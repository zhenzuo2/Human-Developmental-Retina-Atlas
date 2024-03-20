import scanpy as sc
import pandas as pd
import tempfile
import scvi
import scvelo as scv
import os
import sys
x = sys.argv[1]

# Read data
def run_umap_scvi(adata, batch_key="sampleid", labels_key="majorclass"):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata = adata.copy()
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=2000, subset=True
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

df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.csv")
df.index = df['Unnamed: 0'].values

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
cell1 = list(adata[adata.obs.majorclass.isin([x])].obs.index.values)
cell2 = list(df.loc[df.subclass==x,:].index)
cells = list(set(cell1+cell2))
cells = [x for x in cells if x in adata.obs.index]
adata = adata[cells,:]
adata
adata = run_umap_scvi(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/"+x+"_w_NRPC.h5ad"
)
adata.obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/Annotation/"+x+"_w_NRPC.csv"
)