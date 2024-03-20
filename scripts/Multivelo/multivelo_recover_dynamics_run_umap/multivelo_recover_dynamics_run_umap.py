import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import scvi
import tempfile

samples = sys.argv[1]
output_dir = "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_run_umap/"
input_rna_file = "/storage/chentemp/zz4/adult_dev_compare/results/MultiVelo_filtered_cells/adata_rna.h5ad"
input_atac_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_knn_smooth_chrom/"
    + samples
    + ".h5ad"
)
output_file_path = output_dir + samples + ".h5ad"
labels_key = "majorclass"
try:
    os.makedirs(os.path.dirname(output_file_path))
except FileExistsError:
    # directory already exists
    pass

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")
pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 200)
np.set_printoptions(suppress=True)

adata_rna = scv.read(input_rna_file)

adata_atac = scv.read(input_atac_file)

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
if samples == "AC":
    adata = sc.read(
        "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC_w_NRPC.h5ad"
    )
if samples == "BC":
    adata = sc.read(
        "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/BC_w_NRPC.h5ad"
    )
if samples == "Cone":
    adata = sc.read(
        "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/Cone_w_NRPC.h5ad"
    )
if samples == "HC":
    adata = sc.read(
        "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/HC_w_NRPC.h5ad"
    )
if samples == "RGC":
    adata = sc.read(
        "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/RGC_w_NRPC.h5ad"
    )
if samples == "Rod":
    adata = sc.read(
        "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/Rod_w_NRPC.h5ad"
    )
shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata.obs_names))
adata_rna = adata_rna[shared_cells,]

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

adata_rna.obsm["X_umap"] = adata[adata_rna.obs.index].obsm["X_umap"]
adata_rna.write(output_file_path)
