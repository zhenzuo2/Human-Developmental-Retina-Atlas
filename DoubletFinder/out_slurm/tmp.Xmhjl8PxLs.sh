#!/bin/bash
#SBATCH --job-name=tmp.Xmhjl8PxLs
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.Xmhjl8PxLs-%j.out
#SBATCH -e out_slurm/tmp.Xmhjl8PxLs-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript DoubletFinder.R /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object/Multiome_20w1d_FR_nFeature_RNA_1000_3000_MT_5_fitered.rds /storage/singlecell/zz4/fetal_bash/results/DoubletFinder_seurat_object/ /storage/singlecell/zz4/fetal_bash/figures/DoubletFinder_UMAP/ Multiome_20w1d_FR

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
