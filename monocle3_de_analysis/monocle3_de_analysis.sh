input_file="/storage/singlecell/zz4/fetal_snakemake/results/merged_rna/merged_rna.rds"

mkdir /storage/singlecell/zz4/fetal_snakemake/results/AC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_snakemake/results/BC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_snakemake/results/Cone_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_snakemake/results/Rod_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_snakemake/results/HC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_snakemake/results/RGC_monocle3_DE_analysis/
mkdir /storage/singlecell/zz4/fetal_snakemake/results/MG_monocle3_DE_analysis/

meta_file=(
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_AC.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_BC.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_Cone.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_Rod.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_HC.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_RGC.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/RPC.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/PRPC.csv
/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC.csv
)

output_dir=(
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/AC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/BC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Rod_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/HC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/RGC_monocle3_DE_analysis/
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_
/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_
)

for i in "${!meta_file[@]}"
do
	slurmtaco.sh --g01 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Inferior;
	slurmtaco.sh --g01 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Nasal;
	slurmtaco.sh --g01 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Superior;
	slurmtaco.sh --g01 -m 50G -t 1 --30day -- Rscript monocle3_de_analysis.R  $input_file ${meta_file[i]} ${output_dir[i]} All Temporal;
done
