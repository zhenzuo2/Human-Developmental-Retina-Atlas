import scvelo as scv
import anndata
import sys
import os

output_file_path = sys.argv[1]

AC = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/annotated_umap.h5ad"
)
BC = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/annotated_umap.h5ad"
)
Cone = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/annotated_umap.h5ad"
)
Rod = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap.h5ad"
)
HC = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/annotated_umap.h5ad"
)
RGC = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/annotated_umap.h5ad"
)
MG = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap.h5ad"
)

adata = anndata.concat([AC, BC, Cone, Rod, HC, RGC, MG])

isExist = os.path.exists(output_file_path)
if not isExist:
    os.makedirs(output_file_path)

adata.write(output_file_path + "merged_h5ad_adult_annotated.h5ad")
adata.obs.to_csv(output_file_path + "merged_h5ad_adult_annotated_obs.csv")