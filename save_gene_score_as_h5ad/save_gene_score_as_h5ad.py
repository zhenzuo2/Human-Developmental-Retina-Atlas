import scipy
from scipy.io import mmread
import numpy as np
import pandas as pd
import anndata as ad

rownames = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/Gene_score_matrix/gene_score_rownames.txt",
    header=None,
)[0].values

colnames = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/Gene_score_matrix/gene_score_colnames.txt",
    header=None,
)[0].values

gs = scipy.io.mmread(
    "/storage/singlecell/zz4/fetal_bash/results/Gene_score_matrix/gene_score.txt",
)

adata = ad.AnnData(X=gs.transpose().toarray(), obs=colnames, var=rownames)

adata.obs.columns = ["cell_id"]
adata.var.columns = ["gene_name"]

adata.var.index = adata.var.gene_name.values
adata.obs.index = adata.obs.cell_id.values

adata.write(
    "/storage/singlecell/zz4/fetal_bash/results/Gene_score_matrix/gene_score.h5ad"
)