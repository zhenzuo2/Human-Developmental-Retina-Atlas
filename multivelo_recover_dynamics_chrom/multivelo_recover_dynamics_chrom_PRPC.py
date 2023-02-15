import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

pd.options.display.max_columns = None
scv.set_figure_params(dpi=600, dpi_save=600)

adata_rna_file = "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_run_umap_PRPC/adata_umap.h5ad"
adata_atac_file = (
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_knn_smooth_chrom/PRPC.h5ad"
)
output_file = "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.h5ad"

adata_rna = scv.read(adata_rna_file)

All = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_All_feature_selection_FALSE.csv"
)
target = list(set(All.target))
tf = list(set(All.tf))
adata_rna = adata_rna[:, [x for x in adata_rna.var.index if x in target + tf]]

adata_atac = scv.read(adata_atac_file)

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)

print(adata_rna)
print(adata_atac)

sc.pl.umap(adata_rna, color="majorclass")


adata_result = mv.recover_dynamics_chrom(
    adata_rna,
    adata_atac,
    max_iter=5,
    init_mode="invert",
    verbose=False,
    parallel=True,
    save_plot=False,
    rna_only=False,
    fit=True,
    n_anchors=500,
    extra_color_key="majorclass",
    n_jobs=10,
)
adata_result.write(output_file)

################################################
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

adata_result = sc.read(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.h5ad"
)
mv.pie_summary(adata_result)
mv.switch_time_summary(adata_result)

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

mv.velocity_embedding_stream(adata_result, basis="umap", color="Days")
mv.velocity_embedding_stream(adata_result, basis="umap", color="latent_time")
adata_result.obs.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC_obs.csv"
)
adata_result.var.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC_var.csv"
)
