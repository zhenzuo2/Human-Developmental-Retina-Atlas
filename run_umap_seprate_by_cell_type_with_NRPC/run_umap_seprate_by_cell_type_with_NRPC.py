import scanpy as sc
import scvi
import os
import scvelo as scv
import sys
import pandas as pd
import tempfile
import scvi
import os


def run_umap_scvi(adata):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=2000, subset=True
    )
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
    lvae.train(max_epochs=20, n_samples_per_label=100000)
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


input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
label = sys.argv[3]

try:
    os.makedirs(output_file_path)
except FileExistsError:
    # directory already exists
    pass

adata = scv.read("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap/merged_h5ad_adult_annotated_umap_X_scANVI.h5ad")
meta = pd.read_csv(input_file_path)
adata = adata[meta["Unnamed: 0"]]

adata = run_umap_scvi(adata)

adata.obs.to_csv(os.path.join(output_file_path, label + ".csv"))
adata.write(os.path.join(output_file_path, label + ".h5ad"))