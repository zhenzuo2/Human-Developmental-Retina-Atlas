mkdir /storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/
mkdir /storage/chentemp/zz4/adult_dev_compare/figures/subclass_annotation/
export CUDA_VISIBLE_DEVICES=0
adata_query_file="/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult"
cell_type="AC"
adata_ref_file="/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/AC/local.h5ad"
output_result_path="/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/AC_annotation_adult/"
output_fig_path="/storage/chentemp/zz4/adult_dev_compare/figures/subclass_annotation/AC_annotation_adult/"
mkdir ${output_result_path}
mkdir ${output_result_path}
resolution="5";
nfeature="10000"
/storage/chen/home/zz4/miniforge-pypy3/envs/python/bin/python3.9 /storage/chentemp/zz4/adult_dev_compare/scripts/Annotation/subclass_annotation_with_adult/subclass_annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution" "$nfeature" ;

cell_type="BC"
export CUDA_VISIBLE_DEVICES=1
adata_ref_file="/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/BC/local.h5ad"
output_result_path="/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/BC_annotation_adult/"
output_fig_path="/storage/chentemp/zz4/adult_dev_compare/figures/subclass_annotation/BC_annotation_adult/"
CUDA_VISIBLE_DEVICES=1;
resolution="1";
nfeature="10000"
/storage/chen/home/zz4/miniforge-pypy3/envs/python/bin/python3.9 /storage/chentemp/zz4/adult_dev_compare/scripts/Annotation/subclass_annotation_with_adult/subclass_annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution" "$nfeature" ;

cell_type="RGC"
export CUDA_VISIBLE_DEVICES=1
adata_ref_file="/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/RGC/local.h5ad"
output_result_path="/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/RGC_annotation_adult/"
output_fig_path="/storage/chentemp/zz4/adult_dev_compare/figures/subclass_annotation/RGC_annotation_adult/"
CUDA_VISIBLE_DEVICES=2;
resolution="1";
nfeature="10000"
/storage/chen/home/zz4/miniforge-pypy3/envs/python/bin/python3.9 /storage/chentemp/zz4/adult_dev_compare/scripts/Annotation/subclass_annotation_with_adult/subclass_annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution" "$nfeature" ;

cell_type="HC"
export CUDA_VISIBLE_DEVICES=3
adata_ref_file="/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/HC/local.h5ad"
output_result_path="/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/HC_annotation_adult/"
output_fig_path="/storage/chentemp/zz4/adult_dev_compare/figures/subclass_annotation/HC_annotation_adult/"
CUDA_VISIBLE_DEVICES=3;
resolution="1";
nfeature="10000"
/storage/chen/home/zz4/miniforge-pypy3/envs/python/bin/python3.9 /storage/chentemp/zz4/adult_dev_compare/scripts/Annotation/subclass_annotation_with_adult/subclass_annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution" "$nfeature" ;

cell_type="Cone"
export CUDA_VISIBLE_DEVICES=0
adata_ref_file="/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/Photoreceptor/local.h5ad"
output_result_path="/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/Cone_annotation_adult/"
output_fig_path="/storage/chentemp/zz4/adult_dev_compare/figures/subclass_annotation/Cone_annotation_adult/"
CUDA_VISIBLE_DEVICES=3;
resolution="1";
nfeature="10000"
/storage/chen/home/zz4/miniforge-pypy3/envs/python/bin/python3.9 /storage/chentemp/zz4/adult_dev_compare/scripts/Annotation/subclass_annotation_with_adult/subclass_annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution" "$nfeature" ;

cell_type="Rod"
export CUDA_VISIBLE_DEVICES=1
adata_ref_file="/storage/chentemp/zz4/adult_dev_compare/data/adult_reference/Photoreceptor/local.h5ad"
output_result_path="/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation/Rod_annotation_adult/"
output_fig_path="/storage/chentemp/zz4/adult_dev_compare/figures/subclass_annotation/Rod_annotation_adult/"
CUDA_VISIBLE_DEVICES=1;
resolution="1";
nfeature="10000"
/storage/chen/home/zz4/miniforge-pypy3/envs/python/bin/python3.9 /storage/chentemp/zz4/adult_dev_compare/scripts/Annotation/subclass_annotation_with_adult/subclass_annotation_with_adult.py "$adata_query_file" "$cell_type" "$adata_ref_file" "$output_result_path" "$output_fig_path" "$resolution" "$nfeature" ;
