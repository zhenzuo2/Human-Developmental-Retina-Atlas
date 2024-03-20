import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib

matplotlib.rcParams.update({"font.size": 10})
mpl.rcParams["figure.dpi"] = 300

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/Photoreceptor/local.h5ad"
)
adata.X = adata.raw.X
adata.var.index = adata.var.feature_name.values
adata.raw = None
adata_result = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/Cone.h5ad"
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
    'ML_Cone', 'S_Cone'
]

dev_cy = [
    'ML_Cone', 'S_Cone', 'Cone Precursor']

common_cy = [x for x in adult_cy if x in dev_cy]

values = [sc.get.rank_genes_groups_df(adata, group=x).names[0] for x in common_cy]

my_dict = dict(zip(common_cy, values))

sc.pl.dotplot(adata_result, my_dict, "celltype", categories_order=dev_cy)
fig = plt.gcf()
fig.set_size_inches(6, 4)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure3/Cone_dev_dotplot.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)


sc.pl.dotplot(adata[adata.obs.author_cell_type.isin(adult_cy)], my_dict, "author_cell_type", categories_order=adult_cy)
fig = plt.gcf()
fig.set_size_inches(6, 4)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure3/Cone_adult_dotplot.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
