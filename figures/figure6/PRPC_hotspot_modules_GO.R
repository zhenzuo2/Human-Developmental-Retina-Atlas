library(ggplot2)

categories <- 1:33
categories[c(11,22,33)] = ""

categories[1:10] = c("cytoplasmic translation",
"peptide metabolic process",
"translation",
"peptide biosynthetic process",
"amide biosynthetic process",
"amide metabolic process",
"organonitrogen compound biosynthetic process",
"system development ",
"nervous system development",
"multicellular organism development")

categories[12:21] = c("cell cycle process",
"mitotic cell cycle",
"mitotic cell cycle process",
"cell cycle",
"nuclear division",
"chromosome segregation",
"cell division",
"organelle fission",
"regulation of cell cycle process",
"microtubule cytoskeleton organization")

categories[23:32] = c("system development",
"animal organ morphogenesis",
"anatomical structure morphogenesis",
"locomotion",
"oxidative phosphorylation",
"animal organ development",
"extracellular matrix organization",
"extracellular structure organization",
"external encapsulating structure organization",
"aerobic electron transport chain")


categories <- stringr::str_to_title(categories)

values <- 1:30
values[1:10] = c(1.48E-73,
1.19E-24,
2.65E-24,
2.74E-24,
7.35E-23,
7.87E-20,
5.50E-18,
1.48E-15,
1.11E-14,
4.42E-13)
values[12:21] = c(9.87E-44,
2.23E-39,
4.06E-38,
5.64E-37,
4.90E-33,
2.76E-32,
6.49E-31,
1.80E-29,
5.61E-29,
2.17E-27)
values[23:32] = c(1.45477E-06,
1.70173E-06,
4.7257E-05,
6.49195E-05,
9.96825E-05,
0.000111268,
0.000128759,
0.000134901,
0.000148,
0.000164633)


values = -log10(values)
values[c(11,22,33)] <- 0
group = rep(c("Module3","Module2", "Module1"), each = 11)
# Create a data frame
data <- data.frame(Category = categories, Value = values, group = group)
data <- data[33:1,]
# Create the bar plot
ggplot(data, aes(x = Category, y = Value, fill = group)) + geom_bar(stat = "identity") +
  labs(x = "", y = "-log10(Pvalues)") + ggtitle("Biological Process") +
  scale_x_discrete(limits = data$Category) + coord_flip() + theme(text = element_text(size = 30)) +
  scale_fill_manual(values = c("#1f77b4",
                               "#ff7f0e",
                               "#2ca02c",
                               "#d62728",
                               "#9467bd")[3:1],
                    limits = c(#"Module5","Module4", 
                    "Module3", "Module2", "Module1")) + theme(panel.background = element_rect(fill = "transparent"),
                                                                                    plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
                                                                                    panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
                                                                                    legend.box.background = element_rect(fill = "transparent"))+ theme(legend.position = c(0.8, 0.2))

ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure6/PRPC_MG_Modules_GO.svg",
       bg = "transparent", width = 20, height = 13, units = "in")