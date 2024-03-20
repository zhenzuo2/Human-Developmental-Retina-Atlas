names = c("PRPC_1", "PRPC_2", "PRPC_3", "PRPC_4")
for (name in names) {
    BIAS_MODEL = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/bias/models/bias.h5"
    CHROMBPNET_MODEL = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
        name, "/models/chrombpnet.h5", sep = "")
    CHROMBPNET_MODEL_NB = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
        name, "/models/chrombpnet_nobias.h5", sep = "")
    REGIONS = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/",
        name, "/NA_peaks.narrowPeak_no_blacklist.bed", sep = "")
    GENOME = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
    CHROM_SIZES = "/storage/singlecell/zz4/Reference/hg38.chrom.sizes"
    OUT_PREFIX = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_pred_bw/",
        name, "/", sep = "")
    BIGWIG = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
        name, "/auxiliary/data_unstranded.bw", sep = "")

    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_pred_bw/",
        name, ".sh", sep = ""))
    cat(paste("rm -rf ",OUT_PREFIX,sep=""))
    cat("\n")
    cat(paste("mkdir ",OUT_PREFIX,sep=""))
    cat("\n")
    cat(paste("chrombpnet pred_bw \\
  -bm ", BIAS_MODEL, " \\
  -cm ",
        CHROMBPNET_MODEL, " \\
  -cmb ", CHROMBPNET_MODEL_NB, " \\
  -r ",
        REGIONS, " \\
  -g ", GENOME, " \\
  -c ", CHROM_SIZES, " \\
  -op ",
        OUT_PREFIX, " \\
  -bw ", BIGWIG, sep = ""))
    sink()
}