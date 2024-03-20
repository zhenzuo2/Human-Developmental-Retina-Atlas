import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

font = {"family": "Arial", "size": 24}
plt.rc("font", **font)


adata_result = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/AC.h5ad"
)
mv.velocity_embedding_stream(adata_result, basis='umap')
plt.clf()
mv.velocity_embedding_stream(
    adata_result,
    basis="umap",
    color="subclass",
    alpha=1,
    size=20,
    palette={
        "NRPC": "#9467bd",
        "AC Precursor": "#17becf",
        "GABAergic": "#bcbd22",
        "Glycinergic": "#d62728",
        "SACs": "#ff7f0e",
        "dual ACs": "#e377c2",
    },
    title = ""
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/plot_scv.pl.velocity_embedding_stream.tiff",
    transparent=True,
    dpi=300,
)
