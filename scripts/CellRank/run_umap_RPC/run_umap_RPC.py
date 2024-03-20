# Seurat cell cycle
# Human Gene Set: KEGG_CELL_CYCLE https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_CELL_CYCLE
import scanpy as sc
import pandas as pd
import tempfile
import scvi
import scvelo as scv
import os

df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/data/cell_cycle_genes.csv")
cell_cycle_genes = list(df.Gene)

# Read data
def run_umap_scvi(adata, batch_key="sampleid", labels_key="majorclass"):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.obs[labels_key] = adata.obs[labels_key].astype(str)
    adata = adata.copy()
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=2000, subset=True
    )
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, labels_key=labels_key, batch_key=batch_key)
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
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
shared_genes = list(set(adata.var.index) - set(cell_cycle_genes))
adata = adata[:, shared_genes]
adata = adata[adata.obs.majorclass.isin(["NRPC"]),]
adata = run_umap_scvi(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/NRPC.h5ad"
)

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata = adata[adata.obs.majorclass.isin(["PRPC"]),:]
adata = run_umap_scvi(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/PRPC_w_cell_cycle_gene.h5ad"
)

shared_genes = list(set(adata.var.index) - set(cell_cycle_genes))
adata = adata[:, shared_genes]
adata = adata[adata.obs.majorclass.isin(["PRPC"]),:]
adata = run_umap_scvi(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/PRPC.h5ad"
)

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
shared_genes = list(set(adata.var.index) - set(cell_cycle_genes))
adata = adata[:, shared_genes]
adata = adata[adata.obs.majorclass.isin(["PRPC","MG"]),:]
adata = run_umap_scvi(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/PRPC_MG.h5ad"
)


