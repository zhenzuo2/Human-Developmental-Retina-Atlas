import scipy
from scipy.io import mmread
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse

rownames = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score_rownames.txt",
    header=None,
)[0].values

colnames = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score_colnames.txt",
    header=None,
)[0].values

gs = scipy.io.mmread(
    "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score.txt",
)

adata = ad.AnnData(X=gs.transpose().toarray(), obs=colnames, var=rownames)

adata.obs.columns = ["cell_id"]
adata.var.columns = ["gene_name"]

adata.var.index = adata.var.gene_name.values
adata.obs.index = adata.obs.cell_id.values
adata.X = sparse.csr_matrix(adata.X)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/gene_score.h5ad"
)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

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

adata.write("/storage/chentemp/zz4/adult_dev_compare/cellbygene/gene_score.h5ad")
