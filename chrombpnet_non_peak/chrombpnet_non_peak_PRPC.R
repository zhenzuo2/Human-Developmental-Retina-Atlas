sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_non_peaks/chrombpnet_non_peaks_earlyPRPC.sh")
cat("chrombpnet prep nonpeaks -g /storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa -p /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/earlyPRPC/NA_peaks.narrowPeak_no_blacklist.bed -c /storage/singlecell/zz4/Reference/hg38.chrom.sizes -fl /storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json -br /storage/singlecell/zz4/Reference/hg38.blacklist.bed.gz -o /storage/singlecell/zz4/fetal_bash/results/chrombpnet_non_peaks/earlyPRPC/\n")
sink()

sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_non_peaks/chrombpnet_non_peaks_latePRPC.sh")
cat("chrombpnet prep nonpeaks -g /storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa -p /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/latePRPC/NA_peaks.narrowPeak_no_blacklist.bed -c /storage/singlecell/zz4/Reference/hg38.chrom.sizes -fl /storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json -br /storage/singlecell/zz4/Reference/hg38.blacklist.bed.gz -o /storage/singlecell/zz4/fetal_bash/results/chrombpnet_non_peaks/latePRPC/\n")
sink()
