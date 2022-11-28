import scvelo as scv
import glob
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import os

meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/merged_h5ad_adult_annotated_obs.csv"
)
meta.index = meta["Unnamed: 0"].values

adata_atac = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_10w_NR.h5ad"
)
SAMPLES = [
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_14w5d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_10w_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_11w2d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_14w5d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_13W_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_11w2d_FR_2.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_11w2d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_12w3d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_20W2d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_13W_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_19W4d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_23w1d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_12w3d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_19W4d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_20W2d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_14w2d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_16w4d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_16w4d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_23w4d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_20w1d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_23w4d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multi_Fetal_23w1d_FR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_20w1d_NR.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/Multiome_14w2d_NR.h5ad",
]

for f in SAMPLES:
    print("Processing " + f)
    temp = sc.read(f)
    adata_atac = anndata.concat([adata_atac, temp])

adata_rna = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/RNA_velocity/combined.h5ad", cache=False
)
adata_rna.var_names_make_unique()

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_cells = np.intersect1d(shared_cells, meta["Unnamed: 0"])
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

adata_rna.obs["scpred_prediction"] = meta.loc[adata_rna.obs_names, "scpred_prediction"]
adata_rna.obs["Time"] = meta.loc[adata_rna.obs_names, "Time"]
adata_rna.obs["Region"] = meta.loc[adata_rna.obs_names, "Region"]
adata_rna.obs["Days"] = meta.loc[adata_rna.obs_names, "Days"]
adata_rna.obs["subclass"] = meta.loc[adata_rna.obs_names, "subclass"]
adata_rna.obs["majorclass"] = meta.loc[adata_rna.obs_names, "majorclass"]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)

try:
    os.makedirs("/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/")
except FileExistsError:
    # directory already exists
    pass

# Write out filtered cells and prepare to run Seurat WNN --> R script can be found on Github.
adata_rna.obs.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/filtered_cells.csv"
)
adata_rna.write(
    "/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/adata_rna.h5ad"
)
adata_atac.write(
    "/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/adata_atac.h5ad"
)