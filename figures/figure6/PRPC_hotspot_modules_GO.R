library(ggplot2)

categories <- 1:30
categories[c(6,12,18,24,30)] = ""

categories[1:5] = c("system development  ",
"multicellular organismal process",
"multicellular organism development  ",
"nervous system development",
"anatomical structure development ")

categories[7:11] = c("cytoplasmic translation",
"translation",
"peptide biosynthetic process",
"peptide metabolic process",
"amide biosynthetic process")

categories[13:17] = c("cell projection organization",
"multicellular organism development ",
"system development ",
"plasma membrane bounded cell projection organization",
"anatomical structure morphogenesis")

categories[19:23] = c("cell cycle process",
"cell cycle",
"mitotic cell cycle",
"mitotic cell cycle process",
"chromosome segregation")

categories[25:30] = c("system development",
"anatomical structure development",
"multicellular organism development",
"anatomical structure morphogenesis ",
"developmental process")

categories <- stringr::str_to_title(categories)

values <- 1:30
values[1:5] = c(2.11E-23,
2.38E-23,
1.09E-21,
4.38E-20,
7.07E-20)
values[7:11] = c(6.24E-132,
1.20E-81,
3.69E-80,
2.05E-74,
8.87E-74)
values[13:17] = c(5.68E-13,
8.49E-13,
1.20E-12,
1.73E-12,
1.22E-11)
values[19:23] = c(1.32E-69,
2.42E-65,
2.12E-62,
2.37E-57,
1.66E-52)
values[25:30] = c(
1.97E-11,
6.61E-09,
1.12E-08,
3.31E-08,
4.58E-08)

values = -log10(values)
values[c(6,12,18,24,30)] <- 0
group = rep(c("Module5","Module4", "Module3", "Module2", "Module1"), each = 6)
# Create a data frame
data <- data.frame(Category = categories, Value = values, group = group)

# Create the bar plot
ggplot(data, aes(x = Category, y = Value, fill = group)) + geom_bar(stat = "identity") +
  labs(x = "", y = "-log10(Pvalues)") + ggtitle("Biological Process") +
  scale_x_discrete(limits = data$Category) + coord_flip() + theme(text = element_text(size = 30)) +
  scale_fill_manual(values = c("#1f77b4",
                               "#ff7f0e",
                               "#2ca02c",
                               "#d62728",
                               "#9467bd")[5:1],
                    limits = c("Module5","Module4", "Module3", "Module2", "Module1")) + theme(panel.background = element_rect(fill = "transparent"),
                                                                                    plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
                                                                                    panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
                                                                                    legend.box.background = element_rect(fill = "transparent"))+ theme(legend.position = c(0.8, 0.2))

ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure6/PRPC_MG_Modules_GO.svg",
       bg = "transparent", width = 20, height = 13, units = "in")