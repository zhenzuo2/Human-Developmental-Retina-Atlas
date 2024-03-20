import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC_w_NRPC.h5ad")
adata.obs['majorclass'] = adata.obs['majorclass'].astype(str)
adata.obs['subclass'] = adata.obs['subclass'].astype(str)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

markers = {"Module 1":[
"DACH1",
    "EBF1",
    "EBF3",
    "FEZF1",
    "FOXP2",
    "GLI3",
    "HES5",
    "HES6",
    "HEYL",
    "HIF1A",
    "HLF",
    "ISL1",
    "NEUROD1",
    "NEUROD4",
    "NEUROG3",
    "NFIA",
    "NFIB",
    "NFIX",
    "ONECUT1",
    "PRDM13",
    "PROX1",
    "PTF1A",
    "RORA",
    "SATB2",
    "SOX11",
    "SOX2",
    "ST18",
    "TCF4",
    "TSHZ2",
    "ZBTB18",
    "ZNF536"
],"Module 2":[
   "BNC2",
    "CUX2",
    "ESRRG",
    "FOS",
    "FOSB",
    "ID4",
    "JAZF1",
    "LHX2",
    "LHX9",
    "MAF",
    "MEIS2",
    "NR4A2",
    "OTX2",
    "PBX1",
    "PBX3",
    "POU3F3",
    "POU6F2",
    "PRDM8",
    "SOX5",
    "TEAD1",
    "TFAP2C",
    "THRB",
    "TOX",
    "TSHZ1",
    "TSHZ3",
    "ZBTB7C",
    "ZFHX3",
    "ZFPM2",
    "ZNF503"
]}
adata.obs.loc[adata.obs.subclass == "AC Precursor", 'majorclass'] = "AC Precursor"
adata.obs.loc[adata.obs.majorclass == "NRPC", 'majorclass'] = "AC Progenitor"
sc.pl.dotplot(adata,markers,groupby = "majorclass")
fig = plt.gcf()
fig.set_size_inches(20, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure7/plot_AC_GRN_By_Module.tiff",
    dpi=300,
    transparent=True,
    bbox_inches="tight",
)
