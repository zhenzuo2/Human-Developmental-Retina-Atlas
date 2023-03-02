BIAS_MODEL="/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/bias/models/bias.h5"
CHROMBPNET_MODEL="/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/earlyPRPC/models/chrombpnet.h5"
CHROMBPNET_MODEL_NB="/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/earlyPRPC/models/chrombpnet_nobias.h5"
REGIONS="/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/earlyPRPC/NA_peaks.narrowPeak_no_blacklist.bed"
GENOME="/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
CHROM_SIZES="/storage/singlecell/zz4/Reference/hg38.chrom.sizes"
OUT_PREFIX="/storage/singlecell/zz4/fetal_bash/results/chrombpnet_pred_bw/earlyPRPC/"

sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_pred_bw/earlyPRPC.sh")
cat(paste("chrombpnet pred_bw \\
  -bm ", BIAS_MODEL, " \\
  -cm ",CHROMBPNET_MODEL," \\
  -cmb ",
    CHROMBPNET_MODEL_NB, " \\
  -r ", REGIONS, " \\
  -g ", GENOME, " \\
  -c ", CHROM_SIZES, " \\
  -op ", OUT_PREFIX, " \\
", sep = ""))
sink()
