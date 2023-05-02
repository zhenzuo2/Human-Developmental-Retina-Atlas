import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scvi
import scanpy as sc
import os

try:
    os.makedirs("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/")
except FileExistsError:
    # directory already exists
    pass

try:
    os.makedirs("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/")
except FileExistsError:
    # directory already exists
    pass

# load reference
adata_ref = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/data/adult_reference/all/major_clean_major_scvi_Cluster_clean.h5ad"
)
print(adata_ref)
sc.pp.highly_variable_genes(adata_ref, flavor="seurat_v3", n_top_genes=10000,subset = True)
print(adata_ref)

scvi.model.SCVI.setup_anndata(adata_ref, batch_key="sampleid")

arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)

# We train the reference using the standard SCVI workflow, except we add a few non-default parameters that were identified to work well with scArches
vae_ref = scvi.model.SCVI(adata_ref, **arches_params)
vae_ref.train()

adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scVI")
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)

adata_query = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered.h5ad"
)
adata_query = adata_query[:, adata_ref.var_names].copy()
scvi.model.SCVI.prepare_query_anndata(adata_query, vae_ref)

vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category="Unknown",
    labels_key="majorclass",
)

vae_ref_scan.train(max_epochs=100)

adata_ref.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)

vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    vae_ref_scan,
)

vae_q.train(
    max_epochs=250,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
    batch_size = 640
)

adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
adata_query.obs["majorclass"] = vae_q.predict()

adata_full = adata_query.concatenate(adata_ref)
adata_full.obs.batch.cat.rename_categories(["Query", "Reference"], inplace=True)

sc.pp.neighbors(adata_full, use_rep="X_scANVI")
sc.tl.leiden(adata_full)
sc.tl.umap(adata_full)

adata_full.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_wadult_umap_10000.h5ad"
)
adata_query.obs.to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv"
)
