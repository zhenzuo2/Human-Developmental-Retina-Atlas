import matplotlib.pyplot as plt
import matplotlib
import scanpy as sc
import seaborn as sns

matplotlib.rcParams.update({"font.size": 30})
adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata = adata[adata.obs.majorclass.isin(["Rod", "Cone"])]
adata = adata[adata.obs.Region.isin(["Macula","Peripheral"])]
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
#sc.pp.scale(adata,zero_center = False)
adata.obs["Region"] = adata.obs["Region"].replace({"Peripheral":"Periphery"})
adata.obs["Weeks"] = adata.obs.Days.map(
    {
        59: "PCW8",
        70: "PCW10",
        76: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW15",
        103: "PCW15",
        116: "PCW15",
        137: "PCW19",
        141: "PCW19",
        142: "PCW19",
        162: "PCW23",
        165: "PCW23",
    }
)
adata.obs["Weeks"] = (
    adata.obs["Weeks"]
    .astype("category")
    .cat.reorder_categories(["PCW10", "PCW13", "PCW15", "PCW19", "PCW23"])
)

genes = [
    "OTX2",
    "RCVRN",
    "AIPL1",
    "NRL",
    "CRX",
    "PDE6B",
    "NR2E3",
    "ROM1",
    "GNGT2",
    "GNAT1",
    "PDE6H",
]
for gene in genes:
    adata.obs[gene] = adata[:, gene].X.toarray()
    plt.clf()
    sns.violinplot(
        data=adata.obs,
        x="Weeks",
        y=gene,
        hue="Region",
        split=True,
        inner=None,
        palette={"Macula": "#F8766D", "Periphery": "#00BFC4"},
        density_norm='width'
    )
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.xticks(rotation='vertical')
    plt.xlabel('')
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/sns.violinplot_Rod_Cone_"
        + gene
        + ".tiff",
        dpi=600,
        bbox_inches="tight",
        transparent=True,
    )
    plt.clf()

#temp = adata[adata.obs.Weeks.isin(["PCW19","PCW23"])]
#temp
#sc.tl.rank_genes_groups(temp,groupby = "Weeks")
#sc.get.rank_genes_groups_df(temp,group = None).to_csv("/storage/chentemp/zz4/adult_dev_compare/temp/degs.csv")
