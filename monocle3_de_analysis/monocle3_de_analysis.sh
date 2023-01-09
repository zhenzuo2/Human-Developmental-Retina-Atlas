input_file="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"

mkdir /storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/

meta_file=(
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_AC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_BC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Cone.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Rod.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_HC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_RGC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/RPC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/PRPC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/MG.csv
)

output_dir=(
/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/RPC_
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/MG_
)

for i in "${!meta_file[@]}"
do
	slurmtaco.sh --g01 -m 100G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} Early;
	slurmtaco.sh --g01 -m 100G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} Late;
	slurmtaco.sh --g01 -m 100G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All;
done
