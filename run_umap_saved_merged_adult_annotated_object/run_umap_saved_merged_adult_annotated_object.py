import scanpy as sc
import scvi
import os
import scvelo as scv
import sys
import pandas as pd

input_file_path = sys.argv[1]
input_meta_file = sys.argv[2]
output_file_path = sys.argv[3]

try:
   os.makedirs(output_file_path)
except FileExistsError:
   # directory already exists
   pass

adata = scv.read(input_file_path)

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=True)

scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="scpred_prediction")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)

temp = scv.read(input_file_path)

temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
temp.obsm["X_umap"] = adata.obsm["X_umap"]

meta = pd.read_csv(input_meta_file)
meta.index = meta["Unnamed: 0"].values
temp = temp[meta["Unnamed: 0"]]
temp.obs = meta

temp.write(os.path.join(output_file_path, "merged_h5ad_adult_annotated_umap_X_scVI.h5ad"))

scvi.settings.seed = 0
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="scpred_prediction",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100000)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)

temp = scv.read(input_file_path)

temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
temp.obsm["X_umap"] = adata.obsm["X_umap"]
temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]

temp = temp[meta["Unnamed: 0"]]
temp.obs = meta

temp.write(os.path.join(output_file_path, "merged_h5ad_adult_annotated_umap_X_scANVI.h5ad"))