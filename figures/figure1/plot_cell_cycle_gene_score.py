# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
cell_cycle_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", 
"GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", 
"HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", 
"SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", 
"CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", 
"USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8", "HMGB2", 
"CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", 
"CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", 
"SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", 
"ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", 
"CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", 
"CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", 
"ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", 
"CBX5", "CENPA", "CDK2", "CDK4", "CDK6", "CDK7", "CDKN1A", "CDKN1B", 
"STAG1", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "ANAPC10", 
"MAD2L2", "STAG2", "PTTG2", "GADD45G", "DBF4", "YWHAQ", "CHEK1", 
"CHEK2", "CREBBP", "GADD45A", "E2F1", "E2F2", "E2F3", "E2F4", 
"E2F5", "EP300", "ORC6", "ORC3", "CDC26", "ABL1", "ANAPC13", 
"SMC1B", "SFN", "GSK3B", "ANAPC2", "ANAPC4", "HDAC1", "HDAC2", 
"MAD2L1", "SMAD2", "SMAD3", "SMAD4", "MCM3", "MCM7", "MDM2", 
"MYC", "GADD45B", "ATM", "WEE2", "ORC1", "ORC2", "ORC4", "ORC5", 
"FZR1", "ANAPC5", "ANAPC7", "ANAPC11", "PLK1", "ATR", "PRKDC", 
"RAD21", "RB1", "RBL1", "RBL2", "CCND1", "ANAPC1", "SKP1", "SKP2", 
"BUB1B", "TFDP1", "TFDP2", "TGFB1", "TGFB2", "TGFB3", "TP53", 
"SKP1P2", "WEE1", "YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAZ", 
"ZBTB17", "SMC1A", "CDC7", "MAD1L1", "CUL1", "CCNB3", "CDC14B", 
"CDC14A", "CDC23", "CDC16", "CCNA2", "CCNA1", "CCNB1", "CCND2", 
"CCND3", "CCNE1", "CCNH", "PKMYT1", "SMC3", "BUB3", "PTTG1", 
"ESPL1", "CDC25A", "CDC25B", "CDC27", "RBX1"]

sc.tl.score_genes(adata, gene_list=cell_cycle_genes)

scv.pl.umap(
    adata,
    color="score",
    size=50,
    legend_loc="right margin",
    title="",
    vmax = 0.25
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/overall_umap_by_cell_cycle_score.png",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)