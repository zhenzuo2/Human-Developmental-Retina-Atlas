input_file=/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds
meta_file=(
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label/AC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label/BC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label/Cone_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label/Rod_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label/HC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label/RGC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/RPC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/PRPC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/NRPC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/MG_major_sub_class_obs.csv
)
output_dir=(
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label_pseudo_time/
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label_pseudo_time/
/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_with_label_pseudo_time/
/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_with_label_pseudo_time/
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label_pseudo_time/
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label_pseudo_time/
/storage/singlecell/zz4/fetal_bash/results/RPC_annotation_adult_with_label_pseudo_time/
/storage/singlecell/zz4/fetal_bash/results/PRPC_annotation_adult_with_label_pseudo_time/PRPC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/NRPC_annotation_adult_with_label_pseudo_time/NRPC_major_sub_class_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label_pseudo_time/
)

for i in "${!meta_file[@]}"; do
slurmtaco.sh -p gpu -m 100G -t 1 --30day -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript monocle3_pseudotime.R  $input_file ${meta_file[i]} ${output_dir[i]};
done