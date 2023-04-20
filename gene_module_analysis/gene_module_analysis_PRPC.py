import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib

sc.set_figure_params(scanpy=True, dpi=600)
adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad"
)
adata = adata[adata.obs.majorclass == "PRPC"]
Time =pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv"
)
# Time  = pd.read_csv("/Users/zz4/Desktop/PRPC.obs.csv")
adata = adata[Time["Unnamed: 0"]]
adata.obs["latent_time"] = Time["latent_time"].values

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]

cell_cycle_genes = [x.strip() for x in open('/storage/singlecell/zz4/fetal_bash/data/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
adata_cc_genes = adata[:, cell_cycle_genes]
sc.tl.pca(adata_cc_genes)

sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )

adata.layers["counts"] = adata.X
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata)

sc.pp.subsample(adata, n_obs=10000)

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]
adata
adata.obsm["latent_time"] = np.asarray([[x] for x in adata.obs["latent_time"]])

hs = hotspot.Hotspot(
    adata,
    model="danb",
    latent_obsm_key="latent_time",
    umi_counts_obs_key="total_counts",
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations(jobs=10)


hs_genes = hs_results.loc[hs_results.FDR < 0.05].index  # Select genes
local_correlations = hs.compute_local_correlations(
    hs_genes, jobs=10
)  # jobs for parallelization

modules = hs.create_modules(min_gene_threshold=100, core_only=True, fdr_threshold=0.05)
module_scores = hs.calculate_module_scores()


local_correlations.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/PRPC_hs_local_correlations.csv"
)
module_scores.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/PRPC_hs_module_scores.csv"
)
joblib.dump(
    hs,
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/PPRC_hs.create_modules.pkl",
)
joblib.dump(
    adata, "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/PPRC_hs.adata.pkl"
)