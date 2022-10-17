#!/bin/bash
#SBATCH --job-name=tmp.XEVP4qqTUs
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=7-00:00:00
#SBATCH --nodelist=mhgcp-g00
#SBATCH -o out_slurm/tmp.XEVP4qqTUs-%j.out
#SBATCH -e out_slurm/tmp.XEVP4qqTUs-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
velocyto run -o /storage/singlecell/zz4/fetal_bash/results/RNA_velocity/Multi_Fetal_23w1d_FR -@ 12 --samtools-memory 10240 --bcfile /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_23w1d_FR/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -e Multi_Fetal_23w1d_FR /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_23w1d_FR/outs/gex_possorted_bam.bam /storage/singlecell/zz4/Reference/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
