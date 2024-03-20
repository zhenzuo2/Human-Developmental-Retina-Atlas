# Import packages
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import matplotlib

matplotlib.font_manager._load_fontmanager(try_read_cache=False)
plt.rcParams["font.family"] = "Arial"

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC_w_NRPC.h5ad"
)
adata = adata[adata.obs.majorclass =="NRPC"]

sc.pl.umap(adata[adata.obs.Region =="Macula"],color = "Days",frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/plot_AC_NRPC_Days_M.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)

sc.pl.umap(adata[adata.obs.Region =="Peripheral"],color = "Days",frameon=False)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/plot_AC_NRPC_Days_P.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)