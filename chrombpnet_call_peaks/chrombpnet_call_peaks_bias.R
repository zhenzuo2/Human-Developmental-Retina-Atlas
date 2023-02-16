sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks/chrombpnet_call_peaks_bias.sh")
input_bam = "/storage/singlecell/zz4/fetal_bash/data/BPNET_reference/19_D003_Tn5_S1_L001_R1_001.trim.merged.srt.nodup.no_chrM_MT.sort.cleaned.bam"
output_folder = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/bias/"
cat("macs2 callpeak --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits -p 0.01 -t ",
    input_bam, " --outdir ", output_folder, "\n", sep = "")
sink()