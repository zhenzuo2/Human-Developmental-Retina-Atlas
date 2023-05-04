import scanpy as sc
import pandas as pd
import tempfile
import scvi
import scvelo as scv
import os

# Read data
adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw.h5ad")

# Read filtered cells after RNA-seq filtering
meta = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv")

# Find all adult cells
adult_cells = meta.loc[meta.sampleid == "10x3_Lobe_19_D019_NeuN-v7", "Unnamed: 0"]
adult_cells = set(meta["Unnamed: 0"]).intersection(adult_cells)

# Read filtered cells after ATAC-seq filtering
atac_cells = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/results/ATAC_filtered_cells/ATAC_Filtered_cell.csv", sep=" ")
atac_filtered_cells = [string.replace("#", "_") for string in atac_cells.x]

# Find common cells after ATAC-seq and RNA-seq filtering
common_cells = set(meta["Unnamed: 0"]).intersection(atac_filtered_cells)

# Extract common cells and adult cells
adata = adata[list(common_cells) + list(adult_cells)]
meta.index = meta["Unnamed: 0"].values
adata.obs["majorclass"] = meta.loc[adata.obs.index.values, "majorclass"]
adata = adata[
    adata.obs.majorclass.isin(["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]),
]


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

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    temp = scv.read(f)
    os.remove(f)
    temp.obs["leiden"] = adata.obs.leiden
    temp.obsm["X_scVI"] = adata.obsm["X_scVI"]
    temp.obsm["X_umap"] = adata.obsm["X_umap"]
    return temp


adata = run_umap_scvi(adata)
adata.write(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
)
