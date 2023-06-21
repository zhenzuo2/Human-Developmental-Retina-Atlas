library(ggplot2)

categories <- 1:30
categories[c(6, 12, 18, 24,30)] = ""

categories[1:5] = c("system development ",
"nervous system development ",
"multicellular organism development  ",
"multicellular organismal process",
"generation of neurons")

categories[7:11] = c("cell cycle process",
"cell cycle",
"mitotic cell cycle",
"mitotic cell cycle process",
"chromosome segregation")

categories[13:17] = c("system development",
"plasma membrane bounded cell projection organization",
"cell projection organization",
"nervous system development",
"multicellular organism development ")

categories[19:23] = c("cytoplasmic translation",
"translation",
"peptide biosynthetic process",
"peptide metabolic process",
"amide biosynthetic process")

categories[25:29] = c("system development  ",
"anatomical structure development",
"multicellular organism development",
"developmental process",
"anatomical structure morphogenesis")

categories <- stringr::str_to_title(categories)

values <- 1:30
values[1:5] = c(5.08E-12,
1.16E-11,
1.90E-09,
2.42E-09,
1.06E-08)
values[7:11] = c(5.66E-80,
2.38E-74,
4.29E-74,
3.82E-71,
8.47E-69)
values[13:17] = c(4.82E-17,
4.04E-15,
5.72E-15,
3.31E-14,
4.99E-14)
values[19:23] = c(2.16E-141,
9.09E-93,
3.15E-91,
2.79E-85,
1.27E-84)
values[25:29] = c(3.63E-15,
1.35E-12,
4.16E-12,
9.99E-12,
7.76E-11)
values = -log10(values)
values[c(6, 12, 18, 24,30)] <- 0
group = rep(c("Module5", "Module4", "Module3", "Module2", "Module1"), each = 6)
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
    "#9467bd")[5:1],
        limits = c("Module5", "Module4", "Module3", "Module2", "Module1")) + theme(panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"))+ theme(legend.position = c(0.8, 0.2))

ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure6/PRPC_MG_Modules_GO.svg",
    bg = "transparent", width = 20, height = 10, units = "in")