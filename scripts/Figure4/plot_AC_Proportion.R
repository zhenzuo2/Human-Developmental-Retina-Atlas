library(dplyr)
library(ggplot2)
test <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/AC_w_NRPC.csv")
from = c(59, 70, 76, 79, 87, 91, 100, 103, 116, 137, 141, 142, 162, 165)
to = c("PCW8", "PCW10", "PCW10", "PCW10", "PCW13", "PCW13", "PCW15", "PCW15",
    "PCW15", "PCW19", "PCW19", "PCW19", "PCW23", "PCW23")
test$Week <- plyr::mapvalues(test$Days, from = from, to = to)
test$Week <- factor(test$Week, levels = c("PCW8", "PCW10", "PCW13", "PCW15",
    "PCW19", "PCW23"))
test$subclass <- factor(test$subclass, levels = c("dual ACs", "SACs", "Glycinergic",
    "GABAergic", "AC Precursor", "NRPC"))
colordict = c(NRPC = "#9467bd", `AC Precursor` = "#17becf", GABAergic = "#bcbd22",
    Glycinergic = "#d62728", SACs = "#ff7f0e", `dual ACs` = "#e377c2")
result <- test %>%
    group_by(Region, Week, subclass) %>%
    summarise(row_count = n())

write.csv(result,"/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/plot_AC_Proportion.csv")
M = result[result$Region == "Macula", ]
p1 <- ggplot(M, aes(fill = subclass, y = row_count, x = Week)) + geom_bar(position = "fill",
    stat = "identity") + scale_fill_manual(values = colordict[unique(result$subclass)]) +
    scale_y_continuous(labels = scales::percent) + ylab("Percentage") +
    xlab("PCW") + labs(fill = "")

ggsave("/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/plot_AC_Proportion_Macula.tiff",
    plot = p1, width = 6, height = 3, units = "in", dpi = 300) 

P = result[result$Region == "Peripheral", ]
p2 <- ggplot(P, aes(fill = subclass, y = row_count, x = Week)) + geom_bar(position = "fill",
    stat = "identity") + scale_fill_manual(values = colordict[unique(result$subclass)]) +
    scale_y_continuous(labels = scales::percent) + ylab("Percentage") +
    xlab("PCW")

ggsave("/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/plot_AC_Proportion_Peripheral.tiff",
    plot = p2, width = 6, height = 3, units = "in", dpi = 300)
