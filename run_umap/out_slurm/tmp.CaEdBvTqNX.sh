#!/bin/bash
#SBATCH --job-name=tmp.CaEdBvTqNX
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=7-00:00:00
#SBATCH --nodelist=mhgcp-g01
#SBATCH -o out_slurm/tmp.CaEdBvTqNX-%j.out
#SBATCH -e out_slurm/tmp.CaEdBvTqNX-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
python3 run_umap.py /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/ /storage/singlecell/zz4/fetal_bash/results/merged_h5ad/ /storage/singlecell/zz4/fetal_bash/results/scPred_meta/meta.csv

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
