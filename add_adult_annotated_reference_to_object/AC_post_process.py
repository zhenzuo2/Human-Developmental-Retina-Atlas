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
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000.h5ad"
)

scv.pl.umap(adata, color="majorclass", legend_loc="right margin")

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/AC_annotation_adult/AC_major_sub_class.h5ad"
)
scv.pl.umap(adata, color="majorclass", legend_loc="right margin")
scv.pl.umap(adata, color="leiden")


def run_umap_scvi(adata):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, batch_key="sampleid", labels_key="majorclass")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    scvi.settings.seed = 0

    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key="majorclass",
        unlabeled_category="Unknown",
    )
    lvae.train(max_epochs=100, batch_size=640)
    adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp


adata = run_umap_scvi(adata)

scv.pl.umap(adata, color="majorclass", legend_loc="right margin")

scv.pl.umap(adata, color="leiden")
scv.pl.umap(adata, color="sampleid")
adata.obs.loc[~adata.obs.leiden.isin(["12", "24", "8", "34", "35"]), :].to_csv(
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/AC.csv"
)