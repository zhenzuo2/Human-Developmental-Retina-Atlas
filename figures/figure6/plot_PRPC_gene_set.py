import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
adata = scv.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap_PRPC/adata_umap.h5ad")
adata_result = adata_result[adata.obs.index]
adata_result.obsm["X_umap"] = adata.obsm["X_umap"]

sc.pp.highly_variable_genes(
    adata_result, flavor="seurat_v3", n_top_genes=2000, subset=True
)
sc.pp.normalize_total(adata_result)
sc.pp.log1p(adata_result)

genes = [
    "RPL35",
    "RPL10L",
    "RPL22L1",
    "FAU",
    "RPL13A",
    "RPL36",
    "RSL24D1P11",
    "MRPL13",
    "RPSA",
    "RPL10A",
    "RPS27L",
    "RPL26L1",
    "RSL24D1",
    "RPL3",
    "RPL3L",
    "RPL4",
    "RPL5",
    "RPL6",
    "RPL7",
    "RPL7A",
    "RPL8",
    "RPL9",
    "RPL10",
    "RPL11",
    "RPL12",
    "RPL13",
    "RPL15",
    "RPL17",
    "RPL18",
    "RPL18A",
    "RPL19",
    "RPL21",
    "RPL22",
    "RPL23A",
    "RPL24",
    "RPL26",
    "RPL27",
    "RPL30",
    "RPL27A",
    "RPL28",
    "RPL29",
    "RPL31",
    "RPL32",
    "RPL34",
    "RPL35A",
    "RPL36AL",
    "RPL37",
    "RPL37A",
    "RPL38",
    "RPL39",
    "RPL41",
    "RPL36A",
    "RPLP0",
    "RPLP1",
    "RPLP2",
    "RPS2",
    "RPS3",
    "RPS3A",
    "RPS4X",
    "RPS4Y1",
    "RPS5",
    "RPS6",
    "RPS7",
    "RPS8",
    "RPS9",
    "RPS10",
    "RPS11",
    "RPS12",
    "RPS13",
    "RPS15",
    "RPS15A",
    "RPS16",
    "RPS17",
    "RPS18",
    "RPS19",
    "RPS20",
    "RPS21",
    "RPS23",
    "RPS24",
    "RPS25",
    "RPS26",
    "RPS27",
    "RPS27A",
    "RPS28",
    "RPS29",
    "UBA52",
    "RPL14",
    "RPL23",
]
genes = [x for x in genes if x in adata_result.var.index]
print(genes)
sc.tl.score_genes(adata_result, gene_list=genes)
sc.pl.umap(
    adata_result, color="score", vmax=1.5, size=100,title="",frameon = False,colorbar_loc = None
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "Ribosomal Protein Gene Set Enrichment Score.svg",
    bbox_inches="tight",
    backend="cairo",
    transparent=True,
)

genes = ["CDK2","CDK4","CDK6","CDK7","CDKN1A","CDKN1B","STAG1","CDKN1C","CDKN2A","CDKN2B","CDKN2C","CDKN2D","ANAPC10","MAD2L2","STAG2","PTTG2","GADD45G","DBF4","YWHAQ","CHEK1","CHEK2","CREBBP","GADD45A","E2F1","E2F2","E2F3","E2F4","E2F5","EP300","ORC6","ORC3","CDC26","ABL1","ANAPC13","SMC1B","SFN","GSK3B","ANAPC2","ANAPC4","HDAC1","HDAC2","MAD2L1","SMAD2","SMAD3","SMAD4","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MDM2","MYC","GADD45B","ATM","WEE2","ORC1","ORC2","ORC4","ORC5","PCNA","FZR1","ANAPC5","ANAPC7","ANAPC11","PLK1","ATR","PRKDC","RAD21","RB1","RBL1","RBL2","CCND1","ANAPC1","SKP1","SKP2","","","BUB1","BUB1B","TFDP1","TFDP2","TGFB1","TGFB2","TGFB3","TP53","TTK","SKP1P2","","WEE1","YWHAB","YWHAE","YWHAG","YWHAH","YWHAZ","ZBTB17","SMC1A","CDC7","CDC45","MAD1L1","CUL1","CCNB3","CDC14B","CDC14A","CDC23","CDC16","CCNA2","CCNA1","CCNB1","CCND2","CCND3","CCNE1","CCNH","PKMYT1","SMC3","CCNB2","CCNE2","BUB3","PTTG1","ESPL1","CDK1","CDC6","CDC20","CDC25A","CDC25B","CDC25C","CDC27","RBX1"]
genes = [x for x in genes if x in adata_result.var.index]
print(genes)
sc.tl.score_genes(adata_result, gene_list=genes)
sc.pl.umap(
    adata_result, color="score", vmax=0.5, size=100,title="",frameon = False,colorbar_loc = None
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    output_file_path + "KEGG_CELL_CYCLE Gene Set Enrichment Score.svg",
    bbox_inches="tight",
    backend="cairo",
    transparent=True,
)
