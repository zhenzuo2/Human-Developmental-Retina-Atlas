# Download hg38.blacklist.bed.gz from https://github.com/ENCODE-DCC/atac-seq-pipeline and https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks_blacklist/earlyPRPC.sh")
cat("bedtools slop -i /storage/singlecell/zz4/Reference/hg38.blacklist.bed.gz -g /storage/singlecell/zz4/Reference/hg38.chrom.sizes -b 1057 > /storage/singlecell/zz4/Reference/blacklist.bed\n")
cat("bedtools intersect -v -a /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/earlyPRPC/NA_peaks.narrowPeak -b /storage/singlecell/zz4/Reference/blacklist.bed  > /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/earlyPRPC/NA_peaks.narrowPeak_no_blacklist.bed")
sink()

sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks_blacklist/latePRPC.sh")
cat("bedtools slop -i /storage/singlecell/zz4/Reference/hg38.blacklist.bed.gz -g /storage/singlecell/zz4/Reference/hg38.chrom.sizes -b 1057 > /storage/singlecell/zz4/Reference/blacklist.bed\n")
cat("bedtools intersect -v -a /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/latePRPC/NA_peaks.narrowPeak -b /storage/singlecell/zz4/Reference/blacklist.bed  > /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/latePRPC/NA_peaks.narrowPeak_no_blacklist.bed")
sink()