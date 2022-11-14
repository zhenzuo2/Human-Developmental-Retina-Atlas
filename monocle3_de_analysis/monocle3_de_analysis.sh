input_file="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"

meta_file=(
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/PRPC_annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/NRPC_annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap_obs.csv
)

output_dir=(
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_
/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/
)

for i in "${!meta_file[@]}"
do
	slurmtaco.sh --g01 -m 20G -t 1 --30day -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} Early;
	slurmtaco.sh --g01 -m 20G -t 1 --30day -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} Late;
	slurmtaco.sh --g01 -m 20G -t 1 --30day -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All;
done
