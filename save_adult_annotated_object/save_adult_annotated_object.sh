adata_file="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_annotated_umap.h5ad"

cell_type=(
AC
BC
RGC
HC
Cone
RPC
)

meta_cluster_adata_file=(
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult/merged_object.h5ad
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult/merged_object.h5ad
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult/merged_object.h5ad
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult/merged_object.h5ad
/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult/merged_object.h5ad
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult/merged_object.h5ad
)
subclass_reference_file=(
/storage/singlecell/zz4/fetal_bash/data/manual_reference/AC_subclass_annotation.csv
/storage/singlecell/zz4/fetal_bash/data/manual_reference/BC_subclass_annotation.csv
/storage/singlecell/zz4/fetal_bash/data/manual_reference/RGC_subclass_annotation.csv
/storage/singlecell/zz4/fetal_bash/data/manual_reference/HC_subclass_annotation.csv
/storage/singlecell/zz4/fetal_bash/data/manual_reference/Cone_subclass_annotation.csv
/storage/singlecell/zz4/fetal_bash/data/manual_reference/MG_subclass_annotation.csv
)
majorclass_reference_file="/storage/singlecell/zz4/fetal_bash/data/cell_subclass_majorclass_mapping.csv"
output_file_path=(
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult/
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult/
)
output_fig_path=(
/storage/singlecell/zz4/fetal_bash/figures/AC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/figures/BC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/figures/RGC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/figures/HC_annotation_adult/
/storage/singlecell/zz4/fetal_bash/figures/Photoreceptor_annotation_adult/
/storage/singlecell/zz4/fetal_bash/figures/MG_annotation_adult/
)

for i in ${!cell_type[@]}
do
    slurmtaco.sh -p short -m 20G -t 1 -- python3 save_adult_annotated_object.py "$adata_file" "${cell_type[i]}" "${meta_cluster_adata_file[i]}" "${subclass_reference_file[i]}" "$majorclass_reference_file" "${output_file_path[i]}" "${output_fig_path[i]}";
done