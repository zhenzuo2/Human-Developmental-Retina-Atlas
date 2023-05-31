import sys
import scvelo as scv
import scanpy as sc
import anndata
import pandas as pd
import os
import concurrent.futures

def process_item(i):
    # Perform some processing on the item
    # ...
    print(i)
    cell_types = ["AC", "BC", "RGC", "HC", "Photoreceptor", "Photoreceptor"]
    labels = [
        ["amacrine cell"],
        ["OFF-bipolar cell", "ON-bipolar cell"],
        ["retinal ganglion cell"],
        ["retina horizontal cell"],
        ["retinal cone cell"],
        ["retinal rod cell"],
    ]

    cell_type = cell_types[i]
    adata_ref_file = (
        "/storage/singlecell/zz4/fetal_snakemake/data/adult_reference/"
        + cell_type
        + "/local.h5ad"
    )
    if i == 4:
        cell_type = "Cone"
    if i == 5:
        cell_type = "Rod"
    adata_query = sc.read_h5ad(
        "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
    )
    adata_query = adata_query[adata_query.obs.majorclass == cell_type,]
    adata_query.obs.Days = adata_query.obs.Days.astype(float)
    adata_query = adata_query[adata_query.obs.Days > 0]

    adata_ref = scv.read(adata_ref_file)
    adata_ref.X = adata_ref.raw.X
    adata_ref.var.index = list(adata_ref.var.feature_name)
    del adata_query.raw
    del adata_ref.raw
    adata_ref = adata_ref[adata_ref.obs.cell_type.isin(labels[i])]
    adata = anndata.concat([adata_ref, adata_query])

    t = adata.X.toarray()
    pd.DataFrame(data=t, index=adata.obs_names, columns=adata.var_names).to_csv(
        "/storage/singlecell/zz4/fetal_snakemake/results/count_matrx_with_adult/"
        + cell_type
        + ".csv"
    )

num_items = 6

with concurrent.futures.ThreadPoolExecutor() as executor:
    # Submit the process_item function for each item in parallel
    results = executor.map(process_item, range(num_items))

# Access the results
print(list(results))

