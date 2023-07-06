import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({"font.size": 30})
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"
adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC_MG_NRPC/adata_umap.h5ad"
)
adata_result = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
adata_result = adata_result[adata.obs.index]
sc.pp.normalize_total(adata_result)
sc.pp.log1p(adata_result)
adata_result = adata_result[adata.obs.index]
adata_result.obsm["X_umap"] = adata.obsm["X_umap"]
sc.pl.umap(adata_result, color="ATOH7", vmax=1.5, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "ATOH7.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result, color="FOXN4", vmax=2, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "FOXN4.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result, color="NOTCH1", vmax=2, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "NOTCH1.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result, color="HES1", vmax=2, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "HES1.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result, color="HES6", vmax=2, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "HES6.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result, color="HEY2", vmax=1.5, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "HEY2.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "SFRP2",vmax = 3, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "SFRP2.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "DLX1",vmax = 1, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "DLX1.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "DLX2",vmax = 1, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "DLX2.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "ONECUT2",vmax = 1.5, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "ONECUT2.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "NFIA",vmax = 3.5, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "NFIA.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "ASCL1",vmax = 2, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "ASCL1.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "OTX2",vmax = 2, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "OTX2.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "SOX4",vmax = 3, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "SOX4.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "RLBP1",vmax = 2.5, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "RLBP1.png", bbox_inches="tight", transparent=True, dpi=600
)

sc.pl.umap(adata_result,color = "SLC1A3",vmax = 3.5, size=80, frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "SLC1A3.png", bbox_inches="tight", transparent=True, dpi=600
)