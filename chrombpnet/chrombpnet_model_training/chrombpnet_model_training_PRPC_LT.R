names = c("PRPC_1", "PRPC_2", "PRPC_3", "PRPC_4")
for (name in names) {
  dir.create(paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
                 name, "/",sep=""),showWarnings=FALSE)
  bam = paste("/storage/singlecell/zz4/fetal_bash/results/merged_bam/",
              name, "_unsorted_ATAC.sort.cleaned.bam", sep = "")
  g = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
  c = "/storage/singlecell/zz4/Reference/hg38.chrom.sizes"
  p = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/",
            name, "/NA_peaks.narrowPeak_no_blacklist.bed", sep = "")
  n = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_non_peaks/",
            name, "/_negatives.bed", sep = "")
  fl = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json"
  b = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/bias/models/bias.h5"
  output = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
                 name, "/", sep = "")
  sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_model_training/",
             name, ".sh", sep = ""))
  cat(paste("chrombpnet pipeline \\
  -ibam ", bam, " \\
  -d \"ATAC\" \\
  -g ",
            g, " \\
  -c ", c, " \\
  -p ", p, " \\
  -n ", n, " \\
  -fl ",
            fl, " \\
  -b  ", b, " \\
  -o ", output, " \\
", sep = ""))
  sink()
}