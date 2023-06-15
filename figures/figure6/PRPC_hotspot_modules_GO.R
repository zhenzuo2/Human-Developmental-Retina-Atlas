library(ggplot2)

categories <- 1:24
categories[c(6, 12, 18, 24)] = ""

categories[1:5] = c("cytoplasmic translation",
"translation",
"peptide biosynthetic process",
"peptide metabolic process",
"amide biosynthetic process")

categories[7:11] = c("system development ",
"cell projection organization",
"plasma membrane bounded cell projection organization",
"nervous system development ",
"cell morphogenesis")

categories[13:17] = c("mitotic cell cycle",
"cell cycle process",
"mitotic cell cycle process",
"cell cycle",
"chromosome segregation")

categories[19:23] = c("multicellular organismal process",
"system development",
"nervous system development",
"neurogenesis",
"generation of neurons")
categories <- stringr::str_to_title(categories)

values <- 1:23
values[c(6, 12, 18, 24)] <- 0
values[1:5] = c(5.87E-103,
2.39E-52,
5.73E-52,
2.32E-48,
9.89E-47)
values[7:11] = c(6.83E-19,
2.74E-16,
1.08E-15,
1.26E-15,
5.38E-15)
values[13:17] = c(1.15E-70,
1.82E-70,
7.81E-68,
4.98E-65,
1.30E-63)
values[19:23] = c(8.07E-16,
8.39E-15,
1.32E-14,
7.98E-13,
1.28E-12)
values = -log10(values)
values[c(6, 12, 18, 24)] <- 0
group = rep(c("Module1", "Module2", "Module3", "Module4"), each = 6)
# Create a data frame
data <- data.frame(Category = categories, Value = values, group = group)
data = data[24:1, ]
# Create the bar plot
ggplot(data, aes(x = Category, y = Value, fill = group)) + geom_bar(stat = "identity") +
    labs(x = "", y = "-log10(Pvalues)") + ggtitle("Biological Process") +
    scale_x_discrete(limits = data$Category) + coord_flip() + theme(text = element_text(size = 30)) +
    scale_fill_manual(values = c("#636EFA", "#FFA15A", "#00CC96", "#EF553B"),
        limits = c("Module1", "Module2", "Module3", "Module4")) + theme(panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"))+ theme(legend.position = c(0.8, 0.2))

ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure6/PRPC_MG_Modules_GO.svg",
    bg = "transparent", width = 20, height = 10, units = "in")