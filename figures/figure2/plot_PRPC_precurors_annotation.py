import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"
matplotlib.rcParams.update({"font.size": 25})
adata_result = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC_MG_NRPC/adata_umap.h5ad")

adata_result.obs["temp"] = np.nan
adata_result.obs.loc[adata_result.obs.leiden.isin(["14"])&(adata_result.obs.majorclass=="PRPC"), "temp"] = "Neurogenesis"

adata_result.obs.loc[adata_result.obs.leiden.isin(["0","17"])&(adata_result.obs.majorclass=="PRPC"), "temp"] = "MG precursor"

sc.pl.umap(adata_result,color = "temp",size = 50,frameon = False,legend_loc = None,title = "")
fig = plt.gcf()
fig.set_size_inches(10,10)
plt.savefig(
    output_file_path + "RPC_Precursor_annotation.png",
    bbox_inches="tight",
    transparent=True,
    dpi = 600
)

