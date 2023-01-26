adata_query_file="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_annotated_umap.h5ad"
cell_type="AC"
adata_ref_file="/storage/singlecell/zz4/fetal_bash/data/adult_reference/AC/local.h5ad"
output_result_path="/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/AC_annotation_adult/"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/AC_annotation_adult/"
resolution="5";
slurmtaco.sh --g01 -m 20G -t 1 -- python3 annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution";

adata_query_file="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_annotated_umap.h5ad"
cell_type="BC"
adata_ref_file="/storage/singlecell/zz4/fetal_bash/data/adult_reference/BC/local.h5ad"
output_result_path="/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/BC_annotation_adult/"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/BC_annotation_adult/"
CUDA_VISIBLE_DEVICES=1;
resolution="1";
slurmtaco.sh --g01 -m 20G -t 1 -- python3 annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution";

adata_query_file="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_annotated_umap.h5ad"
cell_type="RGC"
adata_ref_file="/storage/singlecell/zz4/fetal_bash/data/adult_reference/RGC/local.h5ad"
output_result_path="/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/RGC_annotation_adult/"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/RGC_annotation_adult/"
CUDA_VISIBLE_DEVICES=2;
resolution="1";
slurmtaco.sh --g01 -m 20G -t 1 -- python3 annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution";

adata_query_file="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_annotated_umap.h5ad"
cell_type="HC"
adata_ref_file="/storage/singlecell/zz4/fetal_bash/data/adult_reference/HC/local.h5ad"
output_result_path="/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/HC_annotation_adult/"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/HC_annotation_adult/"
CUDA_VISIBLE_DEVICES=3;
resolution="1";
slurmtaco.sh --g01 -m 20G -t 1 -- python3 annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution";

cell_type="Cone"
adata_ref_file="/storage/singlecell/zz4/fetal_bash/data/adult_reference/Photoreceptor/local.h5ad"
output_result_path="/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/Cone_annotation_adult/"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/Cone_annotation_adult/"
CUDA_VISIBLE_DEVICES=4;
resolution="1";
slurmtaco.sh --g01 -m 20G -t 1 -- python3 annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution";

cell_type="Rod"
adata_ref_file="/storage/singlecell/zz4/fetal_bash/data/adult_reference/Photoreceptor/local.h5ad"
output_result_path="/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/Rod_annotation_adult/"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/Rod_annotation_adult/"
CUDA_VISIBLE_DEVICES=1;
resolution="1";
slurmtaco.sh --g01 -m 20G -t 1 -- python3 annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution";

cell_type="RPC"
adata_ref_file="/storage/singlecell/zz4/fetal_bash/data/adult_reference/MG/local.h5ad"
output_result_path="/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/MG_annotation_adult/"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/MG_annotation_adult/"
CUDA_VISIBLE_DEVICES=2;
resolution="1";
slurmtaco.sh --g01 -m 20G -t 1 -- python3 annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution";
