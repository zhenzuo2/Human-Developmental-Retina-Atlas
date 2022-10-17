#!/bin/bash
#SBATCH --job-name=tmp.CQQh32gvJN
#SBATCH --partition=short
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.CQQh32gvJN-%j.out
#SBATCH -e out_slurm/tmp.CQQh32gvJN-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript creat_object.R /storage/singlecell/zz4/fetal_snakemake/data/Retina_fetal/ /storage/singlecell/zz4/fetal_bash/results/ArchR/

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
