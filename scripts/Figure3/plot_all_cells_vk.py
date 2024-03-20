import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2
adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata_dynamics.h5ad")

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
vk.plot_projection(color = "majorclass")
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig("/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/plot_all_cells_vk.png",transparent=True,)
