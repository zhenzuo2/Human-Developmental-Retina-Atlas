import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys

samples = sys.argv[1]
adata_atac_file = "/storage/chentemp/zz4/adult_dev_compare/results/MultiVelo_filtered_cells/adata_atac.h5ad"
nn_idx_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_seurat_wnn/"
    + samples
    + "_nn_idx.txt"
)
nn_dist_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_seurat_wnn/"
    + samples
    + "_nn_dist.txt"
)
nn_cells_file = (
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_seurat_wnn/"
    + samples
    + "_nn_cells.txt"
)
output_file_path = (
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_knn_smooth_chrom/"
    + samples
    + ".h5ad"
)

try:
    os.makedirs(os.path.dirname(output_file_path))
except FileExistsError:
    # directory already exists
    pass

adata_atac = scv.read(adata_atac_file)

# Read in Seurat WNN neighbors.
nn_idx = np.loadtxt(
    nn_idx_file,
    delimiter=",",
)
nn_dist = np.loadtxt(
    nn_dist_file,
    delimiter=",",
)
nn_cells = pd.Index(
    pd.read_csv(
        nn_cells_file,
        header=None,
    )[0]
)

# Make sure cell names match.
adata_atac = adata_atac[nn_cells]

print(np.all(nn_cells == adata_atac.obs_names))

mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

adata_atac.write(output_file_path)