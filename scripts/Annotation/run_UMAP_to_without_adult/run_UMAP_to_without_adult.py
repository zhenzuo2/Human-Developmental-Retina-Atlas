import scanpy as sc
import pandas as pd
import tempfile
import scvi
import scvelo as scv
import os

# Read data
adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw.h5ad")


# Read filtered cells after RNA-seq filtering
meta = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/filtered_major_class.csv")

# Find common cells after ATAC-seq and RNA-seq filtering
#common_cells = set(meta["Unnamed: 0"]).intersection(atac_filtered_cells)
common_cells = set(meta["Unnamed: 0"])

# Extract common cells and adult cells
adata = adata[list(common_cells)]
adata.obs_names_make_unique()
adata.var_names_make_unique()
meta.index = meta["Unnamed: 0"].values
meta = meta.loc[list(common_cells),]
adata = adata[list(meta.index)]
adata.obs["majorclass"] = meta.loc[adata.obs.index.values, "majorclass"]
adata = adata[
    adata.obs.majorclass.isin(["AC", "BC", "Cone", "HC", "MG", "RGC", "Rod"]),
]

def run_umap_scvi(adata, batch_key="sampleid", labels_key="majorclass"):
    f = tempfile.mkstemp(suffix=".h5ad")[1]
    adata.obs[labels_key] = adata.obs[labels_key].astype(str)
    adata = adata.copy()
    adata.write(f)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
    )
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, labels_key=labels_key,batch_key = batch_key)
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

adata = run_umap_scvi(adata)
adata.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
)
