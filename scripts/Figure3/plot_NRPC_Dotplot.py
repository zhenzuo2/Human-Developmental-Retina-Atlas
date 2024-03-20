# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib
from PIL import Image
matplotlib.font_manager._load_fontmanager(try_read_cache=False)
plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({'font.size': 10})

adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.csv")
df.index = df['Unnamed: 0'].values

adata = adata[df.index]
df['subclass'] = df['subclass'].replace('NRPC', 'Undetermined')
adata.obs["subclass"] = df['subclass']

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, groupby='subclass', method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=10)
fig = plt.gcf()
fig.set_size_inches(23, 3)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/plot_NRPC_Dotplot.tiff",
    transparent=True,
    bbox_inches="tight",
    dpi = 600
)