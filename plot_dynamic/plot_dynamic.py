import scvelo as scv
import scanpy as sc
import sys
import os

input_path = sys.argv[1]
output_path = sys.argv[2]

try:
   os.makedirs(output_path)
except FileExistsError:
   # directory already exists
   pass

scv.settings.plot_prefix = ""

adata = scv.read(input_path)
adata.obs["scpred_prediction"] = adata.obs["scpred_prediction"].astype(str)
if adata.obs.majorclass.str.contains("MG").any():
    adata.obs.loc[adata.obs.majorclass == "MG", "scpred_prediction"] = "MG"
adata.obs["scpred_prediction"] = adata.obs["scpred_prediction"].astype(
    "category"
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="scpred_prediction",
    title="",
    save=output_path + "scpred_prediction" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata[adata.obs.Region == "Macula"],
    basis="umap",
    color="scpred_prediction",
    title="",
    save=output_path + "Macula_scpred_prediction" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata[adata.obs.Region == "Peripheral"],
    basis="umap",
    color="scpred_prediction",
    title="",
    save=output_path + "Peripheral_scpred_prediction" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="majorclass",
    title="",
    save=output_path + "majorclass" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="subclass",
    title="",
    save=output_path + "subclass" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata, basis="umap", title="", color="Days", save=output_path + "Days" + ".svg"
)
scv.pl.velocity_embedding_stream(
    adata, basis="umap", title="", color="Region", save=output_path + "Region" + ".svg"
)
