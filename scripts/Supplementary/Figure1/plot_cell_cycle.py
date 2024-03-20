import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

adata =sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad")

cy = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/cell_cycle/cell_cycle_meta.csv")

cy.index = cy["Unnamed: 0"].values

adata.obs["CCStage"]=cy.loc[adata.obs.index,"CCStage"]

sc.pl.umap(adata,color = "CCStage",size = 10,frameon = False)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/cy.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
)
