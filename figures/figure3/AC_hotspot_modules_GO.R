library(ggplot2)

categories <- 1:30
categories[c(6, 12, 18, 24,30)] = ""

categories[1:5] = c("cytoplasmic translation",
"translation",
"peptide biosynthetic process",
"amide biosynthetic process",
"peptide metabolic process")

categories[7:11] = c("macromolecule localization",
"cellular localization",
"protein localization",
"cellular macromolecule localization",
"localization")

categories[13:17] = c("regulation of cellular process ",
"central nervous system development",
"nervous system development",
"system development",
"regulation of biological process ")

categories[19:23] = c("organic substance biosynthetic process",
"RNA localization",
"biosynthetic process",
"regulation of nitrogen compound metabolic process",
"positive regulation of nucleobase-containing compound metabolic process")

categories[25:29] = c("regulation of cellular process",
"biological regulation",
"regulation of biological process",
"negative regulation of cellular process",
"regulation of cellular component organization")

categories <- stringr::str_to_title(categories)

values <- 1:30
values[1:5] = c(3.39E-41,
9.38E-24,
1.29E-22,
3.85E-22,
9.05E-21)
values[7:11] = c(7.40E-14,
8.92E-14,
1.84E-12,
2.50E-12,
2.75E-12)
values[13:17] = c(3.23E-17,
1.73E-13,
2.53E-12,
6.80E-12,
2.02E-11)
values[19:23] = c(2.36E-12,
1.16E-11,
1.24E-11,
4.30E-11,
8.08E-11)
values[25:29] = c(1.63E-20,
1.11E-16,
4.95E-16,
1.09E-14,
2.25E-11)
values = -log10(values)
values[c(6, 12, 18, 24,30)] <- 0
group = rep(c("Early", "Early Middle", "Middle", "Late Middle", "Late"), each = 6)
# Create a data frame
data <- data.frame(Category = categories, Value = values, group = group)
data = data[30:1, ]
# Create the bar plot
ggplot(data, aes(x = Category, y = Value, fill = group)) + geom_bar(stat = "identity") +
    labs(x = "", y = "-log10(Pvalues)") + ggtitle("Biological Process") +
    scale_x_discrete(limits = data$Category) + coord_flip() + theme(text = element_text(size = 30)) +
    scale_fill_manual(values = c("#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd"),
        limits = c("Early", "Early Middle", "Middle", "Late Middle", "Late")) + theme(panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"))+ theme(legend.position = c(0.8, 0.2))

ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure3/AC_Modules_GO.svg",
    bg = "transparent", width = 20, height = 10, units = "in")