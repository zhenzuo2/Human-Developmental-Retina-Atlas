import scvelo as scv
import anndata

ldata_file = "/storage/chentemp/zz4/adult_dev_compare/results/RNA_velocity/combined.loom"
output_file = "/storage/chentemp/zz4/adult_dev_compare/results/RNA_velocity/combined.h5ad"

ldata = anndata.read_loom(ldata_file)

ldata.obs.index = [
    x.split(":", 1)[0] + "_" + x.split(":", 1)[1][:-1] + "-1" for x in ldata.obs.index
]
index = list(ldata.obs.index)

ldata.obs.index = index
ldata.raw = ldata

ldata.write(output_file)