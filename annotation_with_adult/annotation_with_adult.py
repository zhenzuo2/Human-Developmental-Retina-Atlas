import sys
import scvelo as scv
import scanpy as sc
import scvi
import anndata
import pandas as pd
import os
import plotly.express as px

# Read both query adata and ref adata
# Ref data was downloaded from https://cellxgene.cziscience.com/collections/af893e86-8e9f-41f1-a474-ef05359b1fb7
adata_query_file = sys.argv[1]
cell_type = sys.argv[2]
adata_ref_file = sys.argv[3]
output_file_path = sys.argv[4]
output_fig_path = sys.argv[5]
resolution = int(sys.argv[6])

try:
   os.makedirs(output_file_path)
except FileExistsError:
   # directory already exists
   pass

try:
   os.makedirs(output_fig_path)
except FileExistsError:
   # directory already exists
   pass

# Read data needed to be annotated
adata_query = scv.read(
    adata_query_file,
    cache=False,
)
adata_query = adata_query[
    adata_query.obs.scpred_prediction == cell_type,
]
adata_query.obs.Days = adata_query.obs.Days.astype(float)

# read reference data
adata_ref = scv.read(
    adata_ref_file
)
adata_ref.X = adata_ref.raw.X
adata_ref.var.index = list(adata_ref.var.feature_name)
del adata_query.raw
del adata_ref.raw

# Run UMAP on query adata to get leiden clusters
batch_key = "batch"
labels_key = ""
adata_query.obs[batch_key] = adata_query.obs[batch_key].astype("str").astype("category")
try:
    sc.pp.highly_variable_genes(
        adata_query,
        flavor="seurat_v3",
        n_top_genes=2000,
        batch_key=batch_key,
        subset=True,
    )

except:
    print("An exception occurred! Fix the batch or run without batch now.")
    sc.pp.highly_variable_genes(
        adata_query, flavor="seurat_v3", n_top_genes=10000, subset=True
    )

scvi.settings.seed = 0
if labels_key == "":
    scvi.model.SCVI.setup_anndata(adata_query, batch_key=batch_key)
else:
    scvi.model.SCVI.setup_anndata(
        adata_query, batch_key=batch_key, labels_key=labels_key
    )
vae = scvi.model.SCVI(adata_query, n_layers=2, n_latent=100, gene_likelihood="nb")
vae.train(use_gpu=True)
adata_query.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata_query, use_rep="X_scVI")
sc.tl.leiden(adata_query, resolution=resolution)
sc.tl.umap(adata_query)

author_cell_type = adata_query.obs.leiden.astype(str).astype("category")
X_scVI = adata_query.obsm["X_scVI"]
X_umap = adata_query.obsm["X_umap"]
leiden = adata_query.obs["leiden"]
cells = adata_query.obs.index
adata_query = scv.read(
    adata_query_file,
    cache=False,
)
adata_query = adata_query[cells, :]
adata_query.obs.Days = adata_query.obs.Days.astype(float)
adata_query.obs["author_cell_type"] = author_cell_type
adata_query.obs["leiden"] = leiden
adata_query.obsm["X_umap"] = X_umap
adata_query.obsm["X_scVI"] = X_scVI
del adata_query.raw

adata_ref.obs["batch"] = adata_ref.obs.sample_uuid
adata_query.obs["sample_source"] = "fetal"
adata_ref.obs["sample_source"] = "adult"

# Concat ref adata and query adata
adata = anndata.concat([adata_ref, adata_query])

batch_key = "batch"
labels_key = ""
adata = adata.copy()
adata.raw = adata  # keep full dimension safe
adata.obs[batch_key] = adata.obs[batch_key].astype("str").astype("category")

# Run UMAP on concated data
try:
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, batch_key=batch_key, subset=True
    )

except:
    print("An exception occurred! Fix the batch or run without batch now.")
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
scvi.settings.seed = 0
if labels_key == "":
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
else:
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, labels_key=labels_key)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=100, gene_likelihood="nb")
vae.train(use_gpu=True)
adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)
sc.tl.umap(adata)

adata.write(
    output_file_path + cell_type + "_merged_object.h5ad"
)

# Plots
df = adata.obs.copy()
df["cell_label"] = adata.obs.sample_source.values
df["x"] = adata.obsm["X_umap"][:, 0]
df["y"] = adata.obsm["X_umap"][:, 1]
fig = px.scatter(
    df,
    x="x",
    y="y",
    color="sample_source",
    width=1200,
    height=1200,
    hover_data=["cell_label"],
    color_discrete_map={"fetal": "#7F7F7F", "adult": "#EF553B"},
)
fig.update_traces(mode="markers", marker_size=5)
fig.update_xaxes(showgrid=False)
fig.update_yaxes(showgrid=False)
fig.update_layout(yaxis_visible=False, yaxis_showticklabels=False)
fig.update_layout(xaxis_visible=False, xaxis_showticklabels=False)
fig.update_layout(
    {
        "plot_bgcolor": "rgba(0, 0, 0, 0)",
        "paper_bgcolor": "rgba(0, 0, 0, 0)",
    }
)
fig.update_layout(legend={"title_text": ""})
fig.update_layout(legend={"itemsizing": "constant"})
fig.update_layout(legend=dict(font=dict(size=10)))
fig.write_html(
    output_fig_path + cell_type + "_subtype_umap_adult_fetal_sample_source.html"
)
fig.write_image(
    output_fig_path + cell_type + "_subtype_umap_adult_fetal_sample_source.svg"
)


df = adata.obs.copy()
df["cell_label"] = adata.obs.author_cell_type.values
df["x"] = adata.obsm["X_umap"][:, 0]
df["y"] = adata.obsm["X_umap"][:, 1]
fig = px.scatter(
    df,
    x="x",
    y="y",
    color="author_cell_type",
    width=1200,
    height=1000,
    hover_data=["cell_label"],
)
fig.update_traces(mode="markers", marker_size=5)
fig.update_xaxes(showgrid=False)
fig.update_yaxes(showgrid=False)
fig.update_layout(yaxis_visible=False, yaxis_showticklabels=False)
fig.update_layout(xaxis_visible=False, xaxis_showticklabels=False)
fig.update_layout(
    {
        "plot_bgcolor": "rgba(0, 0, 0, 0)",
        "paper_bgcolor": "rgba(0, 0, 0, 0)",
    }
)
fig.update_layout(legend={"title_text": ""})
fig.update_layout(legend={"itemsizing": "constant"})
fig.update_layout(legend=dict(font=dict(size=10)))
fig.write_html(
    output_fig_path + cell_type + "_subtype_umap_adult_fetal_author_cell_type.html"
)
fig.write_image(
    output_fig_path + cell_type + "_subtype_umap_adult_fetal_author_cell_type.svg"
)
