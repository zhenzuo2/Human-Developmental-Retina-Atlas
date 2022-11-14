import scvelo as scv
import scanpy as sc
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

scv.settings.plot_prefix = ""

adata = scv.read(input_path)
adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(str)
if adata.obs.majorclass.str.contains("MG").any():
    adata.obs.loc[adata.obs.majorclass == "MG", "scpred_prediction_mode"] = "MG"
adata.obs["scpred_prediction_mode"] = adata.obs["scpred_prediction_mode"].astype(
    "category"
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="scpred_prediction_mode",
    title="",
    save=output_path + "_scpred_prediction_mode" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata[adata.obs.Region == "Macula"],
    basis="umap",
    color="scpred_prediction_mode",
    title="",
    save=output_path + "_Macula_scpred_prediction_mode" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata[adata.obs.Region == "Peripheral"],
    basis="umap",
    color="scpred_prediction_mode",
    title="",
    save=output_path + "Peripheral_scpred_prediction_mode" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="majorclass",
    title="",
    save=output_path + "_majorclass" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="subclass",
    title="",
    save=output_path + "_subclass" + ".svg",
)

scv.pl.velocity_embedding_stream(
    adata, basis="umap", title="", color="Days", save=output_path + "_Days" + ".svg"
)
scv.pl.velocity_embedding_stream(
    adata, basis="umap", title="", color="Region", save=output_path + "_Region" + ".svg"
)
