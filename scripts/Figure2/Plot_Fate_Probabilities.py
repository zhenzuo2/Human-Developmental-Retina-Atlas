import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 8})

output_file_path = "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/"

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata = adata[adata.obs.majorclass.isin(["PRPC"])]

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/PRPC_infer_fate/PRPC_to_terminal_states.csv"
)
df.index = df["Unnamed: 0"].values

adata.obs["MG"] = df.loc[adata.obs.index, "MG"]
adata.obs["NRPC"] = df.loc[adata.obs.index, "NRPC"]

sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)

sc.pl.umap(
    adata,
    color="NRPC",
    vmin=0.99,
    frameon=False,
    title = ""
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_NRPC_fate.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)

sc.pl.umap(
    adata,
    color="MG",
    vmin=0.9,
    frameon=False,
    title = ""
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_MG_fate.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)

adata.obs["Days"] = adata.obs["Days"].astype(float)
sc.pl.umap(
    adata,
    color="Days",
    frameon=False,
    title = ""
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_Days.tiff", bbox_inches="tight", transparent=True, dpi=300
)

sc.pl.umap(
    adata,
    color="leiden",
    frameon=False,
    legend_loc="on data",
    title = ""
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    output_file_path + "PRPC_leiden.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=300,
)
