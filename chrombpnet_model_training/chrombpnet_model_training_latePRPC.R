bam = "/storage/singlecell/zz4/fetal_bash/results/merged_bam/latePRPC_unsorted_ATAC.sort.cleaned.bam"
g = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
c = "/storage/singlecell/zz4/Reference/hg38.chrom.sizes"
p = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/latePRPC/NA_peaks.narrowPeak_no_blacklist.bed"
n = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_non_peaks/latePRPC/_negatives.bed"
fl = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json"
b = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_output/bias/models/bias.h5"
output = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/latePRPC/"
sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_model_training/latePRPC.sh")
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