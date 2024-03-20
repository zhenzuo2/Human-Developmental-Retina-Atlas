import numpy as np
import pandas as pd
import scanpy as sc
import magic
from magic import magic

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score.h5ad")

def smooth(adata):
    magic_op = magic.MAGIC()
    magic_op.set_params(n_jobs=10)
    emt_magic = magic_op.fit_transform(adata.X, genes="all_genes")
    emt_magic = magic_op.transform(genes="all_genes")
    adata.X = emt_magic
    return adata

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata=smooth(adata)
adata.write("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score_imputed.h5ad")

meta = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/Annotation/ALL.csv"
)
meta.index = meta["Unnamed: 0"].values
adata.obs.index = [x.replace("#", "_") for x in adata.obs.index.values]

for x in [
    "sampleid",
    "Time",
    "Region",
    "Days",
    "Data Type",
    "majorclass",
    "leiden",
    "subclass",
    "celltype",
]:
    adata.obs[x] = meta.loc[adata.obs.index.values, x]

adata.write("/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_score_imputed.h5ad")


