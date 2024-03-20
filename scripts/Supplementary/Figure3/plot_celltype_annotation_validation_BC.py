import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib

matplotlib.rcParams.update({"font.size": 20})
mpl.rcParams["figure.dpi"] = 300

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/BC/local.h5ad"
)
adata.X = adata.raw.X
adata.var.index = adata.var.feature_name.values
adata.raw = None
adata_result = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/BC.h5ad"
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
    "BB_GB",
    "DB1",
    "DB2",
    "DB3a",
    "DB3b",
    "DB4a",
    "DB4b",
    "DB5",
    "DB6",
    "FMB",
    "IMB",
    "OFFx",
    "RBC",
]

dev_cy = [
    "BB_GB",
    "DB1",
    "DB2",
    "DB3a",
    "DB3b",
    "DB4a",
    "DB5",
    "DB6",
    "FMB",
    "IMB",
    "OFFx",
    "RBC",
    "BC Precursor",
]

common_cy = [x for x in adult_cy if x in dev_cy]

values = [sc.get.rank_genes_groups_df(adata, group=x).names[0] for x in common_cy]

my_dict = dict(zip(common_cy, values))

sc.pl.dotplot(adata_result, my_dict, "celltype", categories_order=dev_cy)
fig = plt.gcf()
fig.set_size_inches(20, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure3/BC_dev_dotplot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.dotplot(adata, my_dict, "author_cell_type", categories_order=adult_cy)
fig = plt.gcf()
fig.set_size_inches(20, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure3/BC_adult_dotplot.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
