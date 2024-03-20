cys = c("AC_NRPC", "BC_Rod_NRPC", "Cone_NRPC", "HC_NRPC", "RGC_NRPC")
for (cy in cys) {
    dir.create(paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_non_peaks/",cy,sep=""),showWarnings=FALSE)
    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_non_peaks/chrombpnet_non_peaks_",
        cy, ".sh", sep = ""))
    cat(paste("chrombpnet prep nonpeaks \\
  -g /storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa \\
  -p /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/",
        cy, "/NA_peaks.narrowPeak_no_blacklist.bed \\
  -c /storage/singlecell/zz4/Reference/hg38.chrom.sizes \\
  -fl /storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json \\
  -br /storage/singlecell/zz4/Reference/hg38.blacklist.bed.gz \\
  -o /storage/singlecell/zz4/fetal_bash/results/chrombpnet_non_peaks/",
        cy, "/", sep = ""))
    cat("\n")
    sink()
}