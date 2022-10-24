import scvelo as scv
import scanpy as sc
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

scv.settings.plot_prefix = ""

adata = scv.read(input_path)
adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(str)
adata.obs.at[adata.obs.majorclass == "MG", "scpred_prediction_mode"] = "MG"

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="scpred_prediction_mode",
    save=output_path + "_scpred_prediction_mode" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="majorclass",
    save=output_path + "_majorclass" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="subclass",
    save=output_path + "_subclass" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="Days", save=output_path + "_Days" + ".svg"
)
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="Region", save=output_path + "_Region" + ".svg"
)
