sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks_blacklist/bias.sh")
cat("bedtools slop -i /storage/singlecell/zz4/Reference/hg38.blacklist.bed.gz -g /storage/singlecell/zz4/Reference/hg38.chrom.sizes -b 1057 > /storage/singlecell/zz4/Reference/blacklist.bed\n")
cat("bedtools intersect -v -a /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/bias/NA_peaks.narrowPeak -b /storage/singlecell/zz4/Reference/blacklist.bed  > /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/bias/NA_peaks.narrowPeak_no_blacklist.bed")
sink()