bam = "/storage/singlecell/zz4/fetal_bash/data/BPNET_reference/19_D003_Tn5_S1_L001_R1_001.trim.merged.srt.nodup.no_chrM_MT.sort.cleaned.bam"
g = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
c = "/storage/singlecell/zz4/Reference/hg38.chrom.sizes"
p = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/bias/NA_peaks.narrowPeak_no_blacklist.bed"
n = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_non_peaks/bias/_negatives.bed"
fl = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json"
output = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/bias/"
sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_model_training/bias.sh")
cat(paste("chrombpnet bias pipeline \\
  -ibam ", bam, " \\
  -d \"ATAC\" \\
  -g ",
          g, " \\
  -c ", c, " \\
  -p ", p, " \\
  -n ", n, " \\
  -fl ",
          fl, " \\
  -b 0.5 \\
  -o ", output, " \\
", sep = ""))
sink()