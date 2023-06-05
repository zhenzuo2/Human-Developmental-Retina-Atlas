library(ggplot2)
for (celltype in c("Rod", "Cone", "BC", "AC", "RGC", "HC")) {
  time <- read.csv(paste("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/", celltype, "_wadult_monocle3_pseudotime.csv",
                         sep = ""), sep = " ")
  
  time$cellid <- rownames(time)
  
  meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv")
  adult <- time[!(time$cellid %in% meta$X), ]
  time$x <- time$x/mean(adult$x)
  rownames(meta) <- meta$X
  
  meta <- merge(meta, time, by.x = 0, by.y = 0)
  if (celltype == "RGC") {
    meta <- meta[meta$Time %in% c("10w", "11w2d", "12w3d", "13w"),
    ]
  }
  if (celltype == "Rod") {
    meta <- meta[meta$Time %in% c("14w2d", "14w5d", "16w4d", "19w4d",
                                  "20w1d", "20w2d", "23w1d", "23w4d"), ]
  }
  if (celltype == "BC") {
    meta <- meta[meta$Time %in% c("14w5d", "16w4d", "19w4d", "20w1d",
                                  "20w2d", "23w1d", "23w4d"), ]
  }
  ggplot(meta, aes(x = Time, y = x, fill = Region)) + geom_boxplot() +
    xlab("Group") + ylab("Value") + ggtitle("Grouped Boxplot")
  
  mean_values <- aggregate(x ~ Time + Region, data = meta, FUN = mean)
  n1 = table(mean_values$Region)[1]
  n2 = table(mean_values$Region)[2]
  fit1 <- smooth.spline(1:n1, mean_values$x[1:n1], spar = 0.35)
  fit2 <- smooth.spline(1:n2, mean_values$x[(n1 + 1):(n1 + n2)], spar = 0.35)
  mean_values$x2 <- c(predict(fit1)$y, predict(fit2)$y)
  # Plot grouped mean line chart
  ggplot(data = mean_values, aes(x = Time, y = x, group = Region)) +
    geom_line() + geom_point() + labs(x = "Group", y = "Mean Value",
                                      title = "Grouped Mean Line Chart")
  
  p <- ggplot(data = mean_values, aes(x = Time, y = x2, group = Region,
                                      color = Region)) + geom_line(linewidth = 2) + geom_point(size = 5) +
    labs(x = "PCW", y = "Maturation score", title = celltype) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(text = element_text(size = 35))+ theme(rect = element_rect(fill = "transparent"))
  if (!(celltype == "BC")) {
    p = p + theme(legend.position = "none")
  }
  ggsave(paste("/storage/singlecell/zz4/fetal_snakemake/figures/figure1/", celltype, "_pseudotime_by_region_line_chat.svg", sep = ""),
         p, bg = "transparent")
  
}
