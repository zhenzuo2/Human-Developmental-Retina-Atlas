#!/bin/bash
#SBATCH --job-name=tmp.FmmAoIobl1
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH -o out_slurm/tmp.FmmAoIobl1-%j.out
#SBATCH -e out_slurm/tmp.FmmAoIobl1-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
Rscript multivelo_seurat_wnn.R /storage/singlecell/zz4/fetal_bash/results/Pando_merged/seurat_object.rds /storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/RPC.csv RPC /storage/singlecell/zz4/fetal_bash/results/multivelo_seurat_wnn/

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
