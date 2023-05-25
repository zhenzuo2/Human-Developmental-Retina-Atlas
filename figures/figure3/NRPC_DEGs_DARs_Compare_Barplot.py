# Import packages
# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

sc.set_figure_params(transparent=True, fontsize=15)
NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)

NRPC = NRPC[NRPC.obs.leiden!="12"]
NRPC.obs["subclass"] = "Unknown"

clusters = ["2", "13"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["0", "1", "18"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["9", "10"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["4"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["17"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["7", "14"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"

sc.pp.highly_variable_genes(NRPC,flavor='seurat_v3',n_top_genes=5000,subset=True)
sc.pp.normalize_total(NRPC)
sc.pp.log1p(NRPC)

df = NRPC.obs
# define the sample size
sample_size = 6000

# group the dataframe by the 'Group' column
groups = df.groupby("subclass")

# create an empty dataframe to store the sample data
sample_df = pd.DataFrame(columns=df.columns)

for group, data in groups:
    sample_data = data.sample(n=sample_size, replace=True)
    sample_df = pd.concat([sample_df, sample_data])
    
NRPC = NRPC[sample_df.index]
sc.tl.rank_genes_groups(NRPC, "subclass")

gs = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/gene_score.h5ad",
)
gs.obs.index = [x.replace("#", "_") for x in gs.obs.index]
gs = gs[sample_df.index]
gs.obs['subclass'] = NRPC.obs['subclass']
sc.tl.rank_genes_groups(gs, "subclass")

gs_len=[]
g_len =[]
ovp = []
for celltype in ["Cone", "RGC",  "Rod",  "HC","AC", "BC"]:
    print(celltype)
    a = sc.get.rank_genes_groups_df(
        gs, group=celltype, log2fc_min=0.96, pval_cutoff=0.01
    ).names.values
    b = sc.get.rank_genes_groups_df(
        NRPC, group=celltype, log2fc_min=2, pval_cutoff=0.01
    ).names.values
    c = [x for x in a if x in b]
    print(c)
    gs_len = gs_len  + [len(a)]
    g_len = g_len  + [len(b)]
    ovp = ovp  + [len(c)]

def stacked_bar_plot(data):
    set1 = gs_len
    set2 = g_len
    overpaying = ovp
    
    subtracted_array = np.subtract(np.array(set1), np.array(overpaying))
    subtracted = list(subtracted_array)
    x = ["Cone", "RGC",  "Rod",  "HC","AC", "BC"]  # x-axis values
    fig, ax = plt.subplots()
    fig.patch.set_facecolor('none')
    
    ax.bar(x, set1, color='blue', alpha=0.5, label='DAGs')
    ax.bar(x, set2, color='red', alpha=0.5, bottom=subtracted, label='DEGs')
    
    
    ax.set_xlabel('Cell Types')
    ax.set_ylabel('# of Genes')
    ax.set_title('')
    ax.legend()
    ax.grid(False)
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure3/NRPC_DEGs_DARs_Compare_Barplot.svg",
    dpi=600,
    bbox_inches="tight")
    plt.clf()

# Example usage
data = [gs_len, g_len, ovp]
stacked_bar_plot(data)