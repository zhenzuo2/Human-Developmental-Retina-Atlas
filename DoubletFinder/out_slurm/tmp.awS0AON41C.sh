#!/bin/bash
#SBATCH --job-name=tmp.awS0AON41C
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.awS0AON41C-%j.out
#SBATCH -e out_slurm/tmp.awS0AON41C-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript DoubletFinder.R /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object/Multi_Fetal_14w5d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds /storage/singlecell/zz4/fetal_bash/results/DoubletFinder_seurat_object/ /storage/singlecell/zz4/fetal_bash/figures/DoubletFinder_UMAP/ Multi_Fetal_14w5d_FR

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
