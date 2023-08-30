bedtools intersect -a all.bed -b Hu14_17_18_ret_H3K4me2.narrowPeak_hg38_processed.bed H3K27ac.bed -f 0.20 -wa -wb -names H3K4me2 H3K27ac> H3K4me2_H3K27ac_all.bed
bedtools intersect -a p2g.bed -b Hu14_17_18_ret_H3K4me2.narrowPeak_hg38_processed.bed H3K27ac.bed -f 0.20 -wa -wb -names H3K4me2 H3K27ac> H3K4me2_H3K27ac_p2g.bed
library(ggplot2)
# Data
data <- data.frame(Value = c(110856/402224, 155972/402224,94594/402224, 10640/26697, 15002/26697,9655/26697))
data$Group <- c("H3K4me2", "H3K27ac","H3K4me2&H3K27ac","H3K4me2", "H3K27ac","H3K4me2&H3K27ac")
data$Labels <- c("All Peaks", "All Peaks","All Peaks","Linked Peaks", "Linked Peaks","Linked Peaks")

bar_plot <- ggplot(data, aes(x = Group, y = Value, fill = Labels)) + geom_bar(stat = "identity",
    position = "dodge") + labs(title = "", x = "",
    y = "Percentage of Overlapped Peaks") + scale_fill_manual(values = c("blue",
    "red")) + scale_y_continuous(labels = scales::percent) + theme_minimal() +
    theme(text = element_text(size = 25))+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

# Display the plot

# Compare Percentage of Adult Histone Modification Regions
# Histone Post-Translational Modification
# Overlapping with All peaks and Linked Peaks
output_dir <- "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/Save-ALL/"
svg(filename = paste(output_dir, "normalized_effect_barplot.svg", sep = ""),
    width = 10, height = 10)
bar_plot
dev.off()