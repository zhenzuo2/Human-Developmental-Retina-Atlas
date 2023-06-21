import scipy
from scipy.io import mmread
import numpy as np
import pandas as pd
import anndata as ad

rownames = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_rna_PRPC_cell_cycle/cell_cycle_normed_rownames.txt",
    header=None,
)[0].values

colnames = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_rna_PRPC_cell_cycle/cell_cycle_normed_colnames.txt",
    header=None,
)[0].values

input = np.loadtxt("/storage/singlecell/zz4/fetal_snakemake/results/merged_rna_PRPC_cell_cycle/cell_cycle_normed.txt",delimiter=' ')

adata = ad.AnnData(X=input.T, obs=colnames, var=rownames)

adata.obs.columns = ["cell_id"]
adata.var.columns = ["gene_name"]

adata.var.index = adata.var.gene_name.values
adata.obs.index = adata.obs.cell_id.values

adata.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_rna_PRPC_cell_cycle/adata.h5ad"
)