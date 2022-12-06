import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import joblib
from cellrank.tl.kernels import VelocityKernel

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/merged_h5ad_adult_annotated_umap.h5ad"
)
vk = VelocityKernel(adata)
vk.compute_transition_matrix()

adata

joblib.dump(
    vk,
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/merged_h5ad_adult_annotated_umap_vk.pkl",
)