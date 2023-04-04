input_file="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"

mkdir /storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/

meta_file=(
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_AC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_BC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Cone.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Rod.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_HC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_RGC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/RPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/MG.csv
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

region=(
Inferior
Nasal
Superior
Temporal
)

for i in "${!meta_file[@]}"
do
	slurmtaco.sh --g00 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Inferior;
	slurmtaco.sh --g00 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Nasal;
	slurmtaco.sh --g00 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Superior;
	slurmtaco.sh --g00 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Temporal;
done
