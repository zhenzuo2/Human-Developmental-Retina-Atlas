import scanpy as sc
import scvi
import os
import scvelo as scv
import sys
import pandas as pd

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

adata = scv.read(input_file_path)

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=True)

scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="majorclass")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)

scvi.settings.seed = 0
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="majorclass",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=30000)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)

temp = scv.read(input_file_path)

temp.obs["leiden"] = adata.obs.leiden
temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
temp.obsm["X_umap"] = adata.obsm["X_umap"]
temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]

isExist = os.path.exists(output_file_path)
if not isExist:
    os.makedirs(output_file_path)

temp.write(os.path.join(output_file_path, "merged_h5ad_adult_annotated_umap.h5ad"))