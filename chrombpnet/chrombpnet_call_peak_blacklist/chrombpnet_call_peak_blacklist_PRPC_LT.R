# Download hg38.blacklist.bed.gz from https://github.com/ENCODE-DCC/atac-seq-pipeline and https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
names = c("PRPC_1", "PRPC_2", "PRPC_3", "PRPC_4")
for (name in names){
    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks_blacklist/",name,".sh",sep=""))
    cat(paste("bedtools intersect -v -a /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/",name,"/NA_peaks.narrowPeak -b /storage/singlecell/zz4/Reference/blacklist.bed  > /storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/",name,"/NA_peaks.narrowPeak_no_blacklist.bed",sep=""))
    sink()
}