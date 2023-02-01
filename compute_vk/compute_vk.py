import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import joblib
from cellrank.tl.kernels import VelocityKernel

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class_ldata_dynamic.h5ad"
)

adata

vk = VelocityKernel(adata)
vk.compute_transition_matrix()

adata

joblib.dump(
    vk,
    "/storage/singlecell/zz4/fetal_bash/results/vk/merged_raw_filtered_umap_10000_major_sub_class_ldata_dynamic_vk.h5ad",
)