names = c("PRPC_1", "PRPC_2", "PRPC_3", "PRPC_4")
for (name in names) {
    MODEL_H5 = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
        name, "/models/chrombpnet.h5", sep = "")
    REGIONS = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/",
        name, "/NA_peaks.narrowPeak_no_blacklist.bed", sep = "")
    GENOME = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
    CHROM_SIZES = "/storage/singlecell/zz4/Reference/hg38.chrom.sizes"
    OUTPUT_PREFIX = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_contribs_bw/",
        name, "/", sep = "")

    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_contribs_bw/",
        name, ".sh", sep = ""))
    cat(paste("rm -rf ",OUTPUT_PREFIX,sep=""))
    cat("\n")
    cat(paste("mkdir ",OUTPUT_PREFIX,sep=""))
    cat("\n")
    cat(paste("chrombpnet contribs_bw \\
  -m ", MODEL_H5, " \\
  -r ",
        REGIONS, " \\
  -g ", GENOME, " \\
  -c ", CHROM_SIZES, " \\
  -op ",
        OUTPUT_PREFIX, " \\
", sep = ""))
    sink()
}