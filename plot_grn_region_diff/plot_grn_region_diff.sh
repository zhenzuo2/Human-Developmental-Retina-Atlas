mkdir /storage/singlecell/zz4/fetal_bash/figures/Pando_grn/

seurat_pando_object=(
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/Cone_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/Rod_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/HC_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/RGC_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/RPC_Pando/RPC_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/NRPC_seurat_object_Macula.rds
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/Cone_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/Rod_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/HC_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/RGC_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/RPC_Pando/RPC_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_seurat_object_Peripheral.rds
/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/NRPC_seurat_object_Peripheral.rds
)
time_cds=(
/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/RPC_monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/MG_monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/RPC_monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_monocle3.rds
)
labels=(
AC_Macula
BC_Macula
Cone_Macula
Rod_Macula
HC_Macula
RGC_Macula
RPC_Macula
PRPC_Macula
NRPC_Macula
AC_Peripheral
BC_Peripheral
Cone_Peripheral
Rod_Peripheral
HC_Peripheral
RGC_Peripheral
RPC_Peripheral
PRPC_Peripheral
NRPC_Peripheral
)
meta_one=(
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/Cone_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/Rod_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/HC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/RGC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/RPC_Pando/RPC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/NRPC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/Cone_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/Rod_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/HC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/RGC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/RPC_Pando/RPC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/NRPC_modules_meta_Peripheral.csv
)

meta_two=(
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/Cone_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/Rod_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/HC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/RGC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/RPC_Pando/RPC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/NRPC_modules_meta_Peripheral.csv
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/Cone_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/Rod_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/HC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/RGC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/RPC_Pando/RPC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Macula.csv
/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/NRPC_modules_meta_Macula.csv
)

output_dir=(
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/AC_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/Cone_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/Rod_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/HC_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/RGC_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/RPC_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/PRPC_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC_Pando_Macula/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/AC_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/Cone_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/Rod_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/HC_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/RGC_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/RPC_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/PRPC_Pando_Peripheral/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC_Pando_Peripheral/
)

for i in "${!seurat_pando_object[@]}"
do
	slurmtaco.sh -p short -m 50G -t 1 -- Rscript plot_grn_region_diff.R  ${seurat_pando_object[i]} ${time_cds[i]} ${labels[i]} ${meta_one[i]} ${meta_two[i]} ${output_dir[i]};
done


