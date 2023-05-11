import scvelo as scv
import cellrank as cr
import joblib
from cellrank.tl.kernels import VelocityKernel

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata_dynamics.h5ad"
)

adata


vk = VelocityKernel(adata)
vk.compute_transition_matrix()

adata

joblib.dump(
    vk,
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata_dynamics_vk.h5ad"
)