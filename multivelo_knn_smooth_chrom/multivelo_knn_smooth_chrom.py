import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys

adata_atac_file = sys.argv[1]
samples = sys.argv[2]
nn_idx_file = sys.argv[3]
nn_dist_file = sys.argv[4]
nn_cells_file = sys.argv[5]
output_file_path = sys.argv[6]

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