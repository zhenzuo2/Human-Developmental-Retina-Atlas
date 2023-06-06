library(ggplot2)

categories <- 1:24
categories[c(6, 12, 18, 24)] = ""

categories[1:5] = c("cytoplasmic translation", "translation", "peptide biosynthetic process",
    "peptide metabolic process", "amide biosynthetic process")

categories[7:11] = c("cell cycle", "mitotic cell cycle", "mitotic cell cycle process",
    "chromosome segregation", "cell division")

categories[13:17] = c("cell projection organization", "neurogenesis", "plasma membrane bounded cell projection organization",
    "generation of neurons", "neuron differentiation")

categories[19:23] = c("system development", "anatomical structure morphogenesis",
    "nervous system development", "multicellular organism development",
    "multicellular organismal process")
categories <- stringr::str_to_title(categories)

values <- 1:23
values[c(6, 12, 18, 24)] <- 0
values[1:5] = c(2.02e-115, 6.58e-66, 1.74e-64, 2.85e-60, 1.17e-59)
values[7:11] = c(7.51e-80, 1.05e-76, 6.24e-70, 1.72e-62, 2.85e-56)

values[13:17] = c(4.31e-17, 2.3e-16, 2.32e-16, 1.86e-15, 8.85e-15)
values[19:23] = c(1.67e-18, 2.44e-17, 1.29e-15, 2.65e-15, 1.25e-14)
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
    legend.box.background = element_rect(fill = "transparent"))

ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure6/PRPC_MG_Modules_GO.svg",
    bg = "transparent", width = 20, height = 10, units = "in")
