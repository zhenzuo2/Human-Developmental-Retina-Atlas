#!/bin/bash
#SBATCH --job-name=tmp.Vcta9Nuryg
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.Vcta9Nuryg-%j.out
#SBATCH -e out_slurm/tmp.Vcta9Nuryg-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript QC_seurat_apply_filter.R /storage/singlecell/zz4/fetal_bash/results/qc_seurat_object/17W1D_Temporal_retina_min_cell_3_min_features_200.rds /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object/ /storage/singlecell/zz4/fetal_bash/figures/after_qc_seurat_object/ 17W1D_Temporal_retina

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
