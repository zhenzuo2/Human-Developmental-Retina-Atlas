import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib


adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad"
)
adata = adata[adata.obs.majorclass == "PRPC"]
adata.obs.Days = adata.obs.Days.astype(float)
adata = adata[adata.obs.Days > 0]

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]

sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )

adata.layers["counts"] = adata.X

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.tl.pca(adata)

sc.pp.subsample(adata, n_obs=20000)

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]
adata

hs = hotspot.Hotspot(
    adata,
    model="danb",
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_counts",
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations(jobs=12)

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index  # Select genes
local_correlations = hs.compute_local_correlations(
    hs_genes, jobs=12
)  # jobs for parallelization

modules = hs.create_modules(min_gene_threshold=100, core_only=True, fdr_threshold=0.05)
module_scores = hs.calculate_module_scores()

local_correlations.to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/PRPC_hs_local_correlations.csv"
)
module_scores.to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/PRPC_hs_module_scores.csv"
)
joblib.dump(
    hs,
    "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/PPRC_hs.create_modules.pkl",
)
joblib.dump(
    adata, "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/PPRC_hs.adata.pkl"
)