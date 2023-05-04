import scanpy as sc

# Read data
adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
)

# Annotate NRPC
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs.loc[adata.obs.leiden.isin(["8"]), "majorclass"] = "NRPC"
adata.obs.loc[adata.obs.majorclass == "MG", "majorclass"] = "NRPC"

# Annotate MG
adata.obs.loc[adata.obs.leiden.isin(["18", "42"]), "majorclass"] = "MG"

# Annotate PRPC
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
adata.obs.loc[adata.obs.leiden.isin(["45", "4", "2", "53"]), "majorclass"] = "PRPC"

# Remove doublets
adata = adata[~adata.obs.leiden.isin(["56", "39"])]

# Clean cluster 11
cells = adata.obs[adata.obs.leiden == "11"]
remove_cells = cells[~cells.majorclass.isin(["Cone", "NRPC"])].index
adata = adata[~adata.obs.index.isin(remove_cells)]

# Clean cluster 13
cells = adata.obs[adata.obs.leiden == "13"]
remove_cells = cells[~cells.majorclass.isin(["NRPC", "Cone"])].index
adata = adata[~adata.obs.index.isin(remove_cells)]
adata.obs.to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv"
)
adata.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_MG.h5ad"
)
