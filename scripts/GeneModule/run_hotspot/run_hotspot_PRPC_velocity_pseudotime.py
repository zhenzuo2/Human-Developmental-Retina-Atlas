import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib

output_path = "/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/"
adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

adata = adata[adata.obs.majorclass.isin(["PRPC"])]

adata.obs.Days = adata.obs.Days.astype(float)
adata = adata[adata.obs.Days > 0]

meta = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/infer_velocity_pseudotime/PRPC.csv")
meta.index = meta["Unnamed: 0"].values
common_cells = [x for x in meta.index if x in adata.obs.index]
adata = adata[common_cells]
adata.obs["velocity_pseudotime"]=meta.loc[adata.obs.index,"velocity_pseudotime"]

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]

sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
sc.pp.subsample(adata, n_obs=10000, random_state=0)
adata.layers["counts"] = adata.X

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]
adata

adata.obsm["velocity_pseudotime"] = np.asarray([[x] for x in adata.obs["velocity_pseudotime"]])

hs = hotspot.Hotspot(
    adata,
    model="danb",
    layer_key="counts",
    latent_obsm_key="velocity_pseudotime",
    umi_counts_obs_key="total_counts",
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations(jobs=12)

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index  # Select genes
local_correlations = hs.compute_local_correlations(
    hs_genes, jobs=12
)  # jobs for parallelization

modules = hs.create_modules(min_gene_threshold=250, core_only=True, fdr_threshold=0.05)
module_scores = hs.calculate_module_scores()

local_correlations.to_csv(
    output_path + "PRPC_hs_local_correlations_velocity_pseudotime.csv"
)
module_scores.to_csv(
    output_path + "PRPC_hs_module_scores_velocity_pseudotime.csv"
)
joblib.dump(
    hs,
    output_path + "PRPC_hs.create_modules_velocity_pseudotime.pkl",
)
adata.write(output_path + "PRPC_hs.adata_velocity_pseudotime.h5ad"
)
