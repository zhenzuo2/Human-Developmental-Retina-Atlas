import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
adata =sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad")
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
adata.obs["Weeks"] = adata.obs.Days.map(
    {
        70: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW16",
        103: "PCW16",
        116: "PCW16",
        120: "PCW16",
        136: "PCW20",
        137: "PCW20",
        141: "PCW20",
        142: "PCW20",
        162: "PCW23",
        165: "PCW23",
    }
)

adata.obs["Weeks"]  = pd.Categorical(
    adata.obs["Weeks"] , categories=["PCW10", "PCW13", "PCW16","PCW20","PCW23"], ordered=True
)

adata.obs['GPR143'] = [x[0] for x in adata[:,"GPR143"].X.toarray()]
adata.obs['TYR'] = [x[0] for x in adata[:,"TYR"].X.toarray()]
sns.lineplot(data=adata[adata.obs.majorclass.isin(["RGC"])].obs, x="Weeks", y="GPR143", hue="Region")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/GPR143.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()
sns.lineplot(data=adata[adata.obs.majorclass.isin(["BC"])].obs, x="Weeks", y="TYR", hue="Region",)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/TYR.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()