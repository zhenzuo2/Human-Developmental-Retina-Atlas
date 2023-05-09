import scanpy as sc
import multivelo as mv
import matplotlib.pyplot as plt
from glob import glob
import sys
import os

input_path = sys.argv[1]
sample = sys.argv[2]
output_result_path = sys.argv[3]
output_figure_path = sys.argv[4]

try:
   os.makedirs(output_result_path)
except FileExistsError:
   # directory already exists
   pass

try:
   os.makedirs(output_figure_path)
except FileExistsError:
   # directory already exists
   pass

adata_atac = sc.read_10x_mtx(
    input_path + "/outs/filtered_feature_bc_matrix/",
    var_names="gene_symbols",
    cache=False,
    gex_only=False,
)
adata_atac = adata_atac[:, adata_atac.var["feature_types"] == "Peaks"]
adata_atac.obs.index = sample + "_" + adata_atac.obs.index 

print(adata_atac)

if os.path.isfile(input_path + "/outs/peak_annotation.tsv"):
    adata_atac = mv.aggregate_peaks_10x(
        adata_atac,
        input_path + "/outs/peak_annotation.tsv",
        input_path + "/outs/analysis/feature_linkage/feature_linkage.bedpe",
        verbose=True,
    )
if os.path.isfile(input_path + "/outs/atac_peak_annotation.tsv"):
    adata_atac = mv.aggregate_peaks_10x(
        adata_atac,
        input_path + "/outs/atac_peak_annotation.tsv",
        input_path + "/outs/analysis/feature_linkage/feature_linkage.bedpe",
        verbose=True,
    )
plt.hist(adata_atac.X.sum(1), bins=100, range=(0, 100000))
plt.savefig(output_figure_path +sample + ".svg")

print(adata_atac)

adata_atac.write(output_result_path + sample + ".h5ad")