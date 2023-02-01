import sys
import scvelo as scv
import scanpy as sc
import scvi
import anndata
import pandas as pd
import os
import plotly.express as px

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scvi
import scanpy as sc
import tempfile
import scvelo as scv
import os

sc.set_figure_params(dpi=600)

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/Rod_annotation_adult/Rod_major_sub_class.h5ad"
)
scv.pl.umap(adata, color="majorclass", legend_loc="right margin")
scv.pl.umap(adata, color="leiden")

adata.obs["Days"] = adata.obs["Days"].astype(float)
scv.pl.umap(adata, color="Days")

adata = adata[adata.obs.leiden.isin(["1", "2", "10", "11"])]
scv.pl.umap(adata, color="majorclass", legend_loc="right margin")
adata.obs.to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/Rod.csv"
)