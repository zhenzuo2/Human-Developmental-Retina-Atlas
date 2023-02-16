merged_sorted_bam = "/storage/singlecell/zz4/fetal_bash/data/BPNET_reference/19_D003_Tn5_S1_L001_R1_001.trim.merged.srt.nodup.no_chrM_MT.sort.bam"
merged_sorted_cleaned_bam = "/storage/singlecell/zz4/fetal_bash/data/BPNET_reference/19_D003_Tn5_S1_L001_R1_001.trim.merged.srt.nodup.no_chrM_MT.sort.cleaned.bam"
sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_merge_bam/bias.sh")
cat("samtools index ", merged_sorted_bam, "\n", sep = "")
cat("samtools view -o ", merged_sorted_cleaned_bam, " ", merged_sorted_bam,
            " `seq 1 22 | sed 's/^/chr/'`\n", sep = "")
cat("samtools index ", merged_sorted_cleaned_bam, "\n", sep = "")
sink()