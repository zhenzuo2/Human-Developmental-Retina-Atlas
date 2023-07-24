import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt

samples = "PRPC"
pd.options.display.max_columns = None
scv.set_figure_params(dpi=600, dpi_save=600)
input_dir = "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_run_umap"
adata_rna_file = input_dir + "_" + samples + "/adata_umap_full.h5ad"
adata_atac_file = (
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_knn_smooth_chrom/"
    + samples
    + ".h5ad"
)
output_file = (
    "/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/"
    + samples
    + "_all_genes.h5ad"
)

adata_rna = scv.read(adata_rna_file)
genes = pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/figures/figure6/mod_reordered.csv")
genes = list(genes.loc[genes.Module!=-1,"Unnamed: 0"].values)
genes = genes + ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", 
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
genes = list(set(genes))
print(len(genes))
adata_rna = adata_rna[:,adata_rna.var.index.isin(genes)]
adata_atac = scv.read(adata_atac_file)

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna)

print(adata_rna)
print(adata_atac)

sc.pl.umap(adata_rna, color="majorclass")

adata_result = mv.recover_dynamics_chrom(
    adata_rna,
    adata_atac,
    max_iter=5,
    init_mode="invert",
    verbose=True,
    parallel=True,
    save_plot=False,
    rna_only=False,
    fit=True,
    n_anchors=500,
    extra_color_key="majorclass",
    n_jobs=10,
)
adata_result.write(output_file)



