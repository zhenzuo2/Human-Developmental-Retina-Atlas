#!/usr/bin/python
# -*- coding: utf-8 -*-
import scanpy as sc
import scvi
import os
import scvelo as scv
import anndata
import sys
import pandas as pd

input_path = sys.argv[1]
output_path = sys.argv[2]
meta_file = sys.argv[3]

SAMPLES = [
    os.path.join(input_path, folder + "/outs/filtered_feature_bc_matrix.h5")
    for folder in os.listdir(input_path)
]
names = os.listdir(input_path)

# Read data

adata = sc.read_10x_h5(SAMPLES[0])
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.obs.index = [names[0] + "_" + x for x in list(adata.obs.index)]
adata.obs["batch"] = names[0]

names = names[1:]
for (i, f) in enumerate(SAMPLES[1:]):
    print("Processing " + f)
    temp = sc.read_10x_h5(f)
    temp.obs.index = [names[i] + "_" + x for x in list(temp.obs.index)]
    temp.obs["batch"] = names[i]
    temp.var_names_make_unique()
    temp.obs_names_make_unique()
    adata = anndata.concat([adata, temp])

adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.write(os.path.join(output_path, "merged_raw_filtered.h5ad"))

meta = pd.read_csv(meta_file, sep=" ")
adata = adata[meta.index.values]
adata.obs["scpred_prediction"] = meta.loc[adata.obs.index.to_list()].scpred_prediction
adata.obs["Time"] = meta.loc[adata.obs.index.to_list()].Time
adata.obs["Region"] = meta.loc[adata.obs.index.to_list()].Region
adata.obs["Days"] = meta.loc[adata.obs.index.to_list()].Days

adata.obs["scpred_prediction"] = adata.obs["scpred_prediction"].str.replace("MG", "RPC")
adata.write(os.path.join(output_path, "merged_raw_filtered_annotated.h5ad"))

sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", n_top_genes=2000, batch_key="batch", subset=True
)

scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="scpred_prediction")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)

scvi.settings.seed = 0
lvae = scvi.model.SCANVI.from_scvi_model(
    vae, adata=adata, labels_key="scpred_prediction", unlabeled_category="Unknown"
)
lvae.train(max_epochs=20, n_samples_per_label=30000)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)

temp = scv.read(os.path.join(output_path, "merged_raw_filtered_annotated.h5ad"))

temp.obs["leiden"] = adata.obs.leiden
temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
temp.obsm["X_umap"] = adata.obsm["X_umap"]
temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]

temp.write(os.path.join(output_path, "merged_raw_filtered_annotated_umap.h5ad"))
