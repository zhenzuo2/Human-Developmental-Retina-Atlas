import scanpy as sc
import pandas as pd
import tempfile
import scvi
import scvelo as scv
import os
import matplotlib.pyplot as plt
import numpy as np

def run_umap_scvi(adata):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.write(f)
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, batch_key="sampleid", labels_key="majorclass")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp


adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
genes = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/Gene_panel/codebook_0_2024_CP1518.csv")
genes = list(genes.index)
genes = [item for item in genes if not item.startswith("Blank")]

adata = adata[:, genes]
adata = run_umap_scvi(adata)
adata.write("/storage/chentemp/zz4/adult_dev_compare/results/Gene_panel/500.h5ad")
sc.pl.umap(
    adata,
    color="majorclass",
    size=5,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/results/Gene_panel/500.tiff",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)
