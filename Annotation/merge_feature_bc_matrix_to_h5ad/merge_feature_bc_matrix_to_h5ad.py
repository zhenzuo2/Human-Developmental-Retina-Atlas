import scanpy as sc
import os
import scvelo as scv
import anndata
import pandas as pd

input_path = "/storage/singlecell/zz4/fetal_snakemake/data/Retina_fetal/"

try:
   os.makedirs('/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/')
except FileExistsError:
   # directory already exists
   pass

SAMPLES = [
    os.path.join(input_path, folder + "/outs/filtered_feature_bc_matrix.h5")
    for folder in os.listdir(input_path)
]

names = os.listdir(input_path)

# Read data
adata = sc.read_10x_h5(SAMPLES[0])
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.obs.index = [names[0] + "_" + x for x in list(adata.obs.index)]
adata.obs["sampleid"] = names[0]

names = names[1:]
for (i, f) in enumerate(SAMPLES[1:]):
    print("Processing " + f)
    temp = sc.read_10x_h5(f)
    temp.obs.index = [names[i] + "_" + x for x in list(temp.obs.index)]
    temp.obs["sampleid"] = names[i]
    temp.var_names_make_unique()
    temp.obs_names_make_unique()
    adata = anndata.concat([adata, temp])

adata.var_names_make_unique()
adata.obs_names_make_unique()

meta =pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/data/Sample_meta/Retina_fetal_sample_meta.csv")

adata.obs["Time"] = adata.obs.sampleid.map(dict(zip(meta.Samples, meta.Time)))
adata.obs["Region"] = adata.obs.sampleid.map(dict(zip(meta.Samples, meta.Region)))
adata.obs["Days"] = adata.obs.sampleid.map(dict(zip(meta.Samples, meta.Days)))
adata.obs["Data Type"] = adata.obs.sampleid.map(
    dict(zip(meta.Samples, meta["Data Type"]))
)

adata.write("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw.h5ad")
cells = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/DoubletFinder_filtered_cells/DoubletFinder_filtered_cells.csv",
    header=None,
)
adata_filter = adata[
    list(cells[0].values),
]
adata_filter.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered.h5ad"
)