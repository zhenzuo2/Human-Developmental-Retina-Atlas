#!/bin/bash
#SBATCH --job-name=tmp.lmuI1ciI9v
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=7-00:00:00
#SBATCH --nodelist=mhgcp-g00
#SBATCH -o out_slurm/tmp.lmuI1ciI9v-%j.out
#SBATCH -e out_slurm/tmp.lmuI1ciI9v-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
velocyto run -o /storage/singlecell/zz4/fetal_bash/results/RNA_velocity/Multiome_10w_FR -@ 12 --samtools-memory 10240 --bcfile /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_10w_FR/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -e Multiome_10w_FR /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_10w_FR/outs/gex_possorted_bam.bam /storage/singlecell/zz4/Reference/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
