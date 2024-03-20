p2g <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/p2g_unique.bed",
    sep = "\t", header = FALSE)
p2g <- p2g[!duplicated(p2g), ]
p2g_hist <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/compare_histone_p2g.bed",
    sep = "\t", header = FALSE)

p2g_hist_H3K27ac <- p2g_hist[p2g_hist$V4 == "H3K27ac", ]
p2g_hist_H3K27ac <- p2g_hist_H3K27ac[c("V1", "V2", "V3")]
p2g_hist_H3K27ac <- p2g_hist_H3K27ac[!duplicated(p2g_hist_H3K27ac), ]
rownames(p2g_hist_H3K27ac) <- paste(p2g_hist_H3K27ac$V1, p2g_hist_H3K27ac$V2,
    p2g_hist_H3K27ac$V3, sep = "")

p2g_hist_H3K4me2 <- p2g_hist[p2g_hist$V4 == "H3K4me2", ]
p2g_hist_H3K4me2 <- p2g_hist_H3K4me2[c("V1", "V2", "V3")]
p2g_hist_H3K4me2 <- p2g_hist_H3K4me2[!duplicated(p2g_hist_H3K4me2), ]
rownames(p2g_hist_H3K4me2) <- paste(p2g_hist_H3K4me2$V1, p2g_hist_H3K4me2$V2,
    p2g_hist_H3K4me2$V3, sep = "")

common_p2g <- intersect(rownames(p2g_hist_H3K4me2), rownames(p2g_hist_H3K27ac))

all <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/all_peaks.bed",
    sep = "\t", header = FALSE)
all_hist <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/compare_histone_all.bed",
    sep = "\t", header = FALSE)

all_hist_H3K27ac <- all_hist[p2g_hist$V4 == "H3K27ac", ]
all_hist_H3K27ac <- all_hist_H3K27ac[c("V1", "V2", "V3")]
all_hist_H3K27ac <- all_hist_H3K27ac[!duplicated(all_hist_H3K27ac), ]
rownames(all_hist_H3K27ac) <- paste(all_hist_H3K27ac$V1, all_hist_H3K27ac$V2,
    all_hist_H3K27ac$V3, sep = "")

all_hist_H3K4me2 <- all_hist[p2g_hist$V4 == "H3K4me2", ]
all_hist_H3K4me2 <- all_hist_H3K4me2[c("V1", "V2", "V3")]
all_hist_H3K4me2 <- all_hist_H3K4me2[!duplicated(all_hist_H3K4me2), ]
rownames(all_hist_H3K4me2) <- paste(all_hist_H3K4me2$V1, all_hist_H3K4me2$V2,
    all_hist_H3K4me2$V3, sep = "")

common_all <- intersect(rownames(all_hist_H3K4me2), rownames(all_hist_H3K27ac))

value <- c(nrow(all_hist_H3K4me2)/nrow(all), nrow(p2g_hist_H3K4me2)/nrow(p2g),
    nrow(all_hist_H3K27ac)/nrow(all), nrow(p2g_hist_H3K27ac)/nrow(p2g),
    length(common_all)/nrow(all), length(common_p2g)/nrow(p2g))

prop.test(c(nrow(all_hist_H3K4me2), nrow(p2g_hist_H3K4me2)), c(nrow(all),
    nrow(p2g)), p = NULL, alternative = "two.sided", correct = TRUE)

prop.test(c(nrow(all_hist_H3K27ac), nrow(p2g_hist_H3K27ac)), c(nrow(all),
    nrow(p2g)), p = NULL, alternative = "two.sided", correct = TRUE)

prop.test(c(length(common_all), length(common_p2g)), c(nrow(all), nrow(p2g)),
    p = NULL, alternative = "two.sided", correct = TRUE)

# Updated plot_barplot function with group variable
plot_barplot <- function(dataframe, x_col, y_col, group_col, title = "Bar Plot") {
    # Extract data from dataframe
    x_values <- dataframe[[x_col]]
    y_values <- dataframe[[y_col]]
    groups <- dataframe[[group_col]]
    # Create a bar plot with grouped bars
    barplot(matrix(y_values, ncol = length(unique(x_values))), beside = TRUE,
        col = c("#0000FF", "#FF0000"), names.arg = x_values, main = "",
        xlab = "Histone Post-Translational Modification", ylab = "Percentage(%)",
        legend.text = unique(groups), args.legend = list(x = "topright",
            bty = "n"))
}

my_data <- data.frame(
  Category = rep(c("H3K4me2", "H3K27ac", "H3K4me2&H3K27ac"), each = 2),  # Adding a group variable
  Value = value*100,
  Group = rep(c("All OCRs", "Linked OCRs"), 3)  # Adding a group variable
)

# Display the dataframe
print(my_data)
write.csv(my_data,"/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/compare_p2g_adult_histone_modification.csv")

# Testing the plot_barplot function with a group variable
png("/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/compare_p2g_adult_histone_modification.png", width = 1500, height = 1200, bg = "transparent",
    res = 300)
plot_barplot(my_data, "Category", "Value", "Group", "Test Bar Plot with Groups")
dev.off()
