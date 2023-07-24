import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import scvi
import tempfile
samples = "PRPC"

output_dir = "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap"
input_rna_file = "/storage/singlecell/zz4/fetal_snakemake/results/MultiVelo_filtered_cells/adata_rna.h5ad"
input_atac_file = "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_knn_smooth_chrom/"+samples+".h5ad"
output_file_path = output_dir+"_"+samples+"/adata_umap_full.h5ad"
labels_key = "majorclass"
try:
    os.makedirs(os.path.dirname(output_file_path))
except FileExistsError:
    # directory already exists
    pass


def run_umap_scvi(adata,labels_key):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.obs[labels_key]=adata.obs[labels_key].astype(str)
    adata = adata.copy()
    adata.write(f)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=10000, subset=True)
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, labels_key=labels_key)
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
        labels_key=labels_key,
        unlabeled_category="Unknown",
    )
    lvae.train(max_epochs=20)
    adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    
    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp


scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")
pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 200)
np.set_printoptions(suppress=True)

adata_rna = scv.read(input_rna_file)

adata_atac = scv.read(input_atac_file)

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))

len(shared_cells), len(shared_genes)
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

adata_rna = run_umap_scvi(adata_rna,labels_key)
adata_rna.write(output_file_path)