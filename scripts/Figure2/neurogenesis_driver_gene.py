import os
import scvelo as scv
import joblib
import cellrank as cr
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

g = joblib.load(
    "/storage/chentemp/zz4/adult_dev_compare/results/PRPC_infer_fate/PRPC_dynamics_g.pkl"
)

adata = scv.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_ldata.h5ad",
    cache=False,
)
adata = adata[adata.obs.majorclass.isin(["PRPC"])]
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

driver_clusters = ["MG", "NRPC"]
delta_df = g.compute_lineage_drivers(
    lineages=["NRPC"], cluster_key="majorclass", clusters=driver_clusters
)
g.adata.obs["Fate Probabilities to NRPCs"] = g.fate_probabilities["NRPC"].X.flatten()

adata.obs["Fate Probabilities to NRPCs"] = g.adata[adata.obs.index].obs[
    "Fate Probabilities to NRPCs"
]

for x in ["Fate Probabilities to NRPCs"]:
    sc.pl.embedding(
        adata,
        basis="umap",
        color=x,
        color_map="viridis",
        s=10,
        vmin=0.995,
        vmax=1.0,
        frameon=False,
    )
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.title(x, fontsize=20)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/neurogenesis_driver_gene.tiff",
        transparent=True,
        dpi=300,
    )

for x in [
    "LINC01572",
    "MIS18BP1",
    "SMC4",
    "C21orf58",
    "ASPM",
    "APOLD1",
    "TOP2A",
    "CENPE",
    "ECT2",
    "KNL1",
    "ARHGAP11B",
    "KIF14",
    "CEP112",
    "MKI67",
    "TPX2",
    "G2E3",
    "NUSAP1",
    "CKAP5",
    "KIF18A",
    "NUF2",
    "HES6",
    "FOXN4",
]:
    sc.pl.embedding(
        adata,
        basis="umap",
        color=x,
        color_map="viridis",
        s=10,
        vmax="p96",
        frameon=False,
    )
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.title(x, fontsize=40)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/neurogenesis_driver_gene_"
        + x
        + ".tiff",
        transparent=True,
        dpi=300,
    )

for x in ["ATOH7"]:
    sc.pl.embedding(
        adata, basis="umap", color=x, color_map="viridis", s=10, vmax=1, frameon=False
    )
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.title(x, fontsize=40)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/neurogenesis_driver_gene_"
        + x
        + ".tiff",
        transparent=True,
        dpi=300,
    )


def merge_tiffs(input_paths, output_path):
    # Validate input
    if len(input_paths) != 9:
        raise ValueError("Exactly 9 input paths are required.")
    # Open each TIFF image
    images = [Image.open(path) for path in input_paths]
    # Get individual image dimensions
    width, height = images[0].size
    # Create a new image with 3 times the width and 3 times the height
    merged_image = Image.new("RGBA", (width * 3, height * 3))
    # Paste each image into the corresponding position in the new image
    for i in range(3):
        for j in range(3):
            merged_image.paste(images[i * 3 + j], (width * j, height * i))
    # Save the merged image
    merged_image.save(output_path)


genes = ["APOLD1", "ARHGAP11B", "ASPM", "C21orf58", "CENPE", "ECT2", "FOXN4", "KIF14"]
input_paths = [
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/neurogenesis_driver_gene.tiff"
] + [
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/neurogenesis_driver_gene_"
    + x
    + ".tiff"
    for x in genes
]

output_path = "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/neurogenesis_driver_gene_merged.tif"
merge_tiffs(input_paths, output_path)
