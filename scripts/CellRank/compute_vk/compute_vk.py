import scvelo as scv
import cellrank as cr
import joblib


adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata_dynamics.h5ad"
)

adata


vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

adata

joblib.dump(
    vk,
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata_dynamics_vk.h5ad"
)