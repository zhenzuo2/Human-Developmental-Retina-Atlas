#!/bin/bash
#SBATCH --job-name=tmp.gndxEP0QMZ
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.gndxEP0QMZ-%j.out
#SBATCH -e out_slurm/tmp.gndxEP0QMZ-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript QC_seurat.R /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/ Multiome_12w3d_NR 3 200 /storage/singlecell/zz4/fetal_bash/figures/qc_seurat_object/ /storage/singlecell/zz4/fetal_bash/results/qc_seurat_object/

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
