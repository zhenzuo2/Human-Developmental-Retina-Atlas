input_file_path=(
"/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult/AC_major_sub_class.h5ad"
"/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult/BC_major_sub_class.h5ad"
"/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult/Cone_major_sub_class.h5ad"
"/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult/Rod_major_sub_class.h5ad"
"/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult/HC_major_sub_class.h5ad"
"/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult/RGC_major_sub_class.h5ad"
"/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult/RPC_major_sub_class.h5ad"
)

output_file_path=(
"/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/"
"/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/"
"/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/"
"/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/"
"/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/"
"/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/"
"/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/"
)
CUDA_VISIBLE_DEVICES=(
0
1
2
3
0
1
2
)

for i in ${!input_file_path[@]}
do
    CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES[i]"
    slurmtaco.sh --g01 -m 20G -t 1 -- python3 run_umap_saved_adult_annotated_object.py "${input_file_path[i]}" "${output_file_path[i]}";
done