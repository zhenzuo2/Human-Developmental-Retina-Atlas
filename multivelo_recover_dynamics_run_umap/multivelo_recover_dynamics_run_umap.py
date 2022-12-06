import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import scvi

input_rna_file = sys.argv[1]
input_atac_file = sys.argv[2]
output_file_path = sys.argv[3]

try:
    os.makedirs(os.path.dirname(output_file_path))
except FileExistsError:
    # directory already exists
    pass


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

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
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

scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10)

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

adata_rna = run_umap_scvi(adata_rna)

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)

adata_rna.write(output_file_path)