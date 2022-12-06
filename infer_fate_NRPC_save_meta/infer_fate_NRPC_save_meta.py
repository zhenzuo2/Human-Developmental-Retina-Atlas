import scvelo as scv
import scanpy as sc
import pandas as pd

adata = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap/merged_h5ad_adult_annotated_umap_X_scANVI.h5ad"
)
sc.pl.umap(adata, color="scpred_prediction")

NRPC = scv.read("/storage/singlecell/zz4/fetal_bash/results/NRPC_fate.h5ad")
sc.pl.umap(NRPC, color="leiden", legend_loc="on data")


## Annotate AC
sc.pl.umap(NRPC, color="AC")
adata.obs["temp"] = adata.obs.index.isin(
    NRPC[NRPC.obs.leiden.isin(["2", "3", "17", "19"])].obs.index
).astype(str)
sc.pl.umap(adata, color="temp")
pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(["2", "3", "17", "19"])].obs.index].obs,
        adata[adata.obs.scpred_prediction == "AC"].obs,
    ]
).to_csv("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_AC.csv")

## Annotate Cone
sc.pl.umap(NRPC, color="Cone")
adata.obs["temp"] = adata.obs.index.isin(
    NRPC[NRPC.obs.leiden.isin(["5"])].obs.index
).astype(str)
sc.pl.umap(adata, color="temp")
pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(["5"])].obs.index].obs,
        adata[adata.obs.scpred_prediction == "Cone"].obs,
    ]
).to_csv("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Cone.csv")

## Annnotate BC
sc.pl.umap(NRPC, color="BC")
adata.obs["temp"] = adata.obs.index.isin(
    NRPC[NRPC.obs.leiden.isin(["6", "7", "9", "15", "22"])].obs.index
).astype(str)
sc.pl.umap(adata, color="temp")
pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(["6", "7", "9", "15", "22"])].obs.index].obs,
        adata[adata.obs.scpred_prediction == "BC"].obs,
    ]
).to_csv("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_BC.csv")

## Annotate Rod
sc.pl.umap(NRPC, color="Rod")
adata.obs["temp"] = adata.obs.index.isin(
    NRPC[NRPC.obs.leiden.isin(["21", "1", "11"])].obs.index
).astype(str)
sc.pl.umap(adata, color="temp")
pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(["21", "1", "11"])].obs.index].obs,
        adata[adata.obs.scpred_prediction == "Rod"].obs,
    ]
).to_csv("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Rod.csv")

## Annotate RGC
sc.pl.umap(NRPC, color="RGC")
adata.obs["temp"] = adata.obs.index.isin(
    NRPC[NRPC.obs.leiden.isin(["16"])].obs.index
).astype(str)
sc.pl.umap(adata, color="temp")
pd.concat(
    [
        adata[NRPC[NRPC.obs.leiden.isin(["16"])].obs.index].obs,
        adata[adata.obs.scpred_prediction == "RGC"].obs,
    ]
).to_csv("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_RGC.csv")

## Annotate HC
pd.concat(
    [
        adata[adata.obs.scpred_prediction == "HC"].obs,
    ]
).to_csv("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_HC.csv")