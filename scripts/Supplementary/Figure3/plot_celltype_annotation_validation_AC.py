import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib

matplotlib.rcParams.update({"font.size": 20})
mpl.rcParams["figure.dpi"] = 300

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/AC/local.h5ad"
)
adata.X = adata.raw.X
adata.var.index = adata.var.feature_name.values
adata.raw = None
adata_result = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC.h5ad"
)

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000, subset=True)

common_genes = list(set(adata.var.index) & set(adata_result.var.index))

adata = adata[:, common_genes]
adata_result = adata_result[:, common_genes]

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.normalize_total(adata_result)
sc.pp.log1p(adata_result)

sc.tl.rank_genes_groups(adata, "author_cell_type")

adult_cy = [
    "AC0",
    "AC1",
    "AC2",
    "AC3",
    "AC4",
    "AC5",
    "AC6",
    "AC7",
    "AC8",
    "AC9",
    "AC10",
    "AC11",
    "AC12",
    "AC13",
    "AC14",
    "AC15",
    "AC16",
    "AC17",
    "AC18",
    "AC19",
    "AC20",
    "AC21",
    "AC22",
    "AC23",
    "AC24",
    "AC25",
    "AC26",
    "AC27",
    "AC28",
    "AC29",
    "AC30",
    "AC31",
    "AC32",
]

dev_cy = [
    "AC0",
    "AC1",
    "AC2",
    "AC3",
    "AC4",
    "AC5",
    "AC6",
    "AC7",
    "AC8",
    "AC9",
    "AC10",
    "AC11",
    "AC12",
    "AC13",
    "AC14",
    "AC15",
    "AC16",
    "AC17",
    "AC18",
    "AC19",
    "AC20",
    "AC21",
    "AC22",
    "AC23",
    "AC24",
    "AC25",
    "AC26",
    "AC28",
    "AC29",
    "AC30",
    "AC31",
    "AC Precursor",
]

common_cy = [x for x in adult_cy if x in dev_cy]

values = [sc.get.rank_genes_groups_df(adata, group=x).names[0] for x in common_cy]

my_dict = dict(zip(common_cy, values))

sc.pl.dotplot(adata_result, my_dict, "celltype", categories_order=dev_cy)
fig = plt.gcf()
fig.set_size_inches(20, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure3/AC_dev_dotplot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.dotplot(adata, my_dict, "author_cell_type", categories_order=adult_cy)
fig = plt.gcf()
fig.set_size_inches(20, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure3/AC_adult_dotplot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
