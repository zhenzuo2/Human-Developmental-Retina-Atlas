input_path=(
/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap_ldata.h5ad
)
output_path=(
/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_umap_ldata_dynamic.h5ad
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/annotated_umap_ldata_dynamic.h5ad
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/annotated_umap_ldata_dynamic.h5ad
/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/annotated_umap_ldata_dynamic.h5ad
/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap_ldata_dynamic.h5ad
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/annotated_umap_ldata_dynamic.h5ad
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/annotated_umap_ldata_dynamic.h5ad
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap_ldata_dynamic.h5ad
)

for i in ${!input_path[@]}
do
    slurmtaco.sh -p gpu -m 20G -t 1 -- python3 recover_dynamic.py "${input_path[i]}" "${output_path[i]}";
done;
