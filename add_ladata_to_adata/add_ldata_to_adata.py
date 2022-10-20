import sys
import scvelo as scv
import anndata

adata_file = sys.argv[1]
ldata_file = sys.argv[2]
output_file = sys.argv[3]

adata = scv.read(adata_file)
ldata = anndata.read_loom(ldata_file)

ldata.obs.index = [
    x.split(":", 1)[0] + "_" + x.split(":", 1)[1][:-1] + "-1" for x in ldata.obs.index
]
index = list(ldata.obs.index)

ldata.obs.index = index
ldata.raw = ldata

adata = scv.utils.merge(adata, ldata)
adata.write(output_file)