import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/"

adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/subclass_annotation/AC_annotation_adult/AC_merged_object.h5ad"
)

adata.obs["temp"] = np.nan
adata.obs.loc[adata.obs.sample_source == "adult", "temp"] = adata.obs.loc[
    adata.obs.sample_source == "adult", "author_cell_type"
]
sc.pl.umap(
    adata,
    color="temp",
    legend_loc="on data",
    frameon=False,
    title="",
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "AC_annotation.svg",
    bbox_inches="tight",
)
