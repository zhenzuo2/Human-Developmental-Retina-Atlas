import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
import joblib
import tempfile
import scvi


def run_umap_scvi(adata, labels_key):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.obs[labels_key] = adata.obs[labels_key].astype(str)
    adata = adata.copy()
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, labels_key=labels_key)
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    scvi.settings.seed = 0
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key=labels_key,
        unlabeled_category="Unknown",
    )
    lvae.train(max_epochs=20)
    adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_scANVI"] = adata.obsm["X_scANVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp


import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
import joblib

adata = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG_ldata.h5ad"
)

meta = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_AC.csv"
)
adata = adata[meta["Unnamed: 0"]]

meta_file = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/AC_subtype.csv"
)
meta_file.index = meta_file["Unnamed: 0"]
common_cells = [x for x in meta_file.index if x in adata.obs.index]
meta_file = meta_file.loc[common_cells, :]
adata.obs["subclass"] = adata.obs["majorclass"].astype(str)
adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)

for sublcass in set(meta_file.subclass):
    adata.obs.loc[
        meta_file.loc[meta_file.subclass == sublcass, :].index, "subclass"
    ] = sublcass


for majorclass in set(meta_file.majorclass):
    adata.obs.loc[
        meta_file.loc[meta_file.majorclass == majorclass, :].index, "majorclass"
    ] = majorclass

adata.obs["Days"] = adata.obs["Days"].astype(float)
adata = run_umap_scvi(adata, "majorclass")
adata.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/cellrank_result/AC_adata.h5ad",
)
sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", n_top_genes=2000, subset=True
)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(
    adata, basis="umap", legend_fontsize=12, title="", color="majorclass"
)

vk = VelocityKernel(adata)
vk.compute_transition_matrix()

ck = ConnectivityKernel(adata).compute_transition_matrix()

combined_kernel = 0.8 * vk + 0.2 * ck

g = GPCCA(combined_kernel)
print(g)
g.compute_schur(n_components=10)
g.plot_spectrum()


terminal_states = (
    adata.obs["majorclass"].replace("NRPC", np.nan).replace("AC Precursor", np.nan)
)
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities()

temp = pd.DataFrame(g.adata.obsm["to_terminal_states"])
temp.columns = list(g.adata.obsm["to_terminal_states"].names)
temp.index = g.adata.obs.index

temp.to_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cellrank_result/AC_to_terminal_states.csv"
)

joblib.dump(
    g,
    "/storage/singlecell/zz4/fetal_snakemake/results/cellrank_result/AC_g.h5ad",
)
