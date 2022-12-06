import scanpy as sc
import scvi
import os
import scvelo as scv
import sys
import pandas as pd

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

try:
   os.makedirs(output_file_path)
except FileExistsError:
   # directory already exists
   pass

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

temp = scv.read(input_file_path)

temp.obs["leiden"] = adata.obs.leiden
temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
temp.obsm["X_umap"] = adata.obsm["X_umap"]


temp.obs.to_csv(os.path.join(output_file_path, "annotated_umap_obs.csv"))
temp.write(os.path.join(output_file_path, "annotated_umap.h5ad"))