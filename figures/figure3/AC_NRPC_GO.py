import scanpy as sc
import pandas as pd
meta =pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_AC.csv")
meta =meta.loc[meta.majorclass =="NRPC", ]
adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad")
adata = adata[meta["Unnamed: 0"]]
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=True)
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

sc.pp.neighbors(adata,use_rep='X')
sc.tl.leiden(adata)
sc.tl.umap(adata)
sc.pl.umap(adata,color = "leiden",legend_loc = "on data")
sc.tl.rank_genes_groups(adata, 'leiden')
df = sc.get.rank_genes_groups_df(adata,group="2")
df.loc[(df.pvals_adj<0.01)&(df.logfoldchanges>0),"names"]
for gene in df.loc[(df.pvals_adj<0.01)&(df.logfoldchanges>0),"names"]:
    print(gene)