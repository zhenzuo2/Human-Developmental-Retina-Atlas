import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/AC.h5ad"
)
meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/AC_NRPC_SCENT.csv"
)
meta.index = meta.X.values
adata.obs["time"] = meta.loc[adata.obs.index, "ccat"]

meta_file = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/AC_subtype.csv"
)
meta_file.index = meta_file["Unnamed: 0"]
common_cells = [x for x in meta_file.index if x in adata.obs.index]
meta_file = meta_file.loc[common_cells, :]
adata.obs["subclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)

for sublcass in set(meta_file.subclass):
    adata.obs.loc[
        meta_file.loc[meta_file.subclass == sublcass, :].index, "subclass"
    ] = sublcass
for majorclass in set(meta_file.majorclass):
    adata.obs.loc[
        meta_file.loc[meta_file.majorclass == majorclass, :].index, "majorclass"
    ] = majorclass

sc.pl.umap(
    adata,
    color="time",
    cmap="plasma",
    vmax=0.11,
    vmin=0.05,
    frameon=False,
    title="",
    size=20,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_ccat.svg",
    dpi=600,
)
plt.clf()

sns.boxplot(x="Days", y="time", data=adata.obs[adata.obs.Region == "Macula"])
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_ccat_Macula_boxplot_Days.svg",
    dpi=600,
)
plt.clf()

sns.boxplot(x="Days", y="time", data=adata.obs[adata.obs.Region == "Peripheral"])
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_ccat_Peripheral_boxplot_Days.svg",
    dpi=600,
)
plt.clf()

sns.boxplot(
    x="majorclass",
    y="time",
    data=adata.obs,
    order=["NRPC", "AC Precursor", "SACs", "GABAergic", "Glycinergic", "dual ACs"],
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_ccat_boxplot_majorclass.svg",
    dpi=600,
)
plt.clf()
