import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib

output_path = "/storage/singlecell/zz4/fetal_snakemake/results/hotspot/PRPC/monocle3/"
adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad"
)

adata = adata[adata.obs.majorclass.isin(["PRPC","MG"])]
adata.obs.Days = adata.obs.Days.astype(float)
adata = adata[adata.obs.Days > 0]

meta = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_MG_monocle3_pseudotime.csv",sep = " ")
adata.obs["monocle3"]=meta.loc[adata.obs.index,"x"]

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]

sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )

adata.layers["counts"] = adata.X

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]
adata

adata.obsm["monocle3"] = np.asarray([[x] for x in adata.obs["monocle3"]])

hs = hotspot.Hotspot(
    adata,
    model="danb",
    layer_key="counts",
    latent_obsm_key="monocle3",
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
    output_path + "PRPC_hs_local_correlations.csv"
)
module_scores.to_csv(
    output_path + "PRPC_hs_module_scores.csv"
)
joblib.dump(
    hs,
    output_path + "PPRC_hs.create_modules.pkl",
)
joblib.dump(
    adata, output_path + "PPRC_hs.adata.pkl"
)