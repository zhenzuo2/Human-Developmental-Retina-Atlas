import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"
matplotlib.rcParams.update({"font.size": 25})
adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC.h5ad")
adata.raw = None
adata_result = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC_MG_NRPC/adata_umap.h5ad")

adata_result.obs["temp"] = np.nan
adata_result.obs.loc[
    adata[adata.obs.leiden.isin(["10"])].obs.index, "temp"
] = "Neurogenesis"

adata_result.obs.loc[adata[adata.obs.leiden.isin(["0"])].obs.index, "temp"] = "MG precursor"

sc.pl.umap(adata_result,color = "temp",size = 50,frameon = False,legend_loc = None)
fig = plt.gcf()
fig.set_size_inches(10,10)
plt.savefig(
    output_file_path + "RPC_Precursor_annotation.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

