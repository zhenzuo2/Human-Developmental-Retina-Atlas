import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
adata =sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad")
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

markers = ["DHRS3","RBP4","RDH10","FGF8","STRA6","CYP26A1","CYP26B1","CYP26C1","RXRG","RXRA","RXRB","RARA","RARB","RARG","ALDH1A1","ALDH1A2","ALDH1A3"]
for gene in markers:
    vmax = np.quantile(adata[:,gene].X.toarray(),q = 0.999)
    sc.pl.umap(adata[adata.obs.Region == "Macula"],color = gene,vmax = vmax,size = 20)
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/"+gene+"_Macula.png",
        dpi=600,
        bbox_inches="tight",
    )
    plt.clf()
    sc.pl.umap(adata[adata.obs.Region == "Peripheral"],color = gene,vmax = vmax,size = 20)
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/"+gene+"_Peripheral.png",
        dpi=600,
        bbox_inches="tight",
    )
    plt.clf()

sc.pl.violin(adata[adata.obs.majorclass == "PRPC"],keys ="ALDH1A1",groupby="Region")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/ALDH1A1_violin.png",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

adata.obs["majorclass"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC/MG",
        "PRPC": "RPC/MG",
        "MG": "RPC/MG",
    }
)
adata.obs["majorclass_Region"] = adata.obs["majorclass"].astype(str)+" "+adata.obs["Region"].astype(str)

sc.pl.violin(adata,keys ="HMX1",groupby="majorclass_Region",rotation = 90)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/HMX1_violin.png",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.violin(adata[adata.obs.majorclass == "PRPC"],keys ="CYP26A1",groupby="Region")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/CYP26A1_violin.png",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.violin(adata[adata.obs.majorclass == "PRPC"],keys ="FGF8",groupby="Region")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/FGF8_violin.png",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

sc.pl.violin(adata[adata.obs.majorclass == "Rod"],keys ="RXRA",groupby="Region")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/RXRA_violin.png",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()


adata.obs["cell_type_Region"] = [m+"_"+n for m,n in zip(adata.obs["majorclass"],adata.obs["Region"])]
sc.pl.dotplot(adata, markers, groupby='cell_type_Region', dendrogram=False)
fig = plt.gcf()
fig.set_size_inches(8, 8)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure5/dotplot_Macula_Peripheral.svg",
    dpi=600,
    bbox_inches="tight",
)
plt.clf()

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