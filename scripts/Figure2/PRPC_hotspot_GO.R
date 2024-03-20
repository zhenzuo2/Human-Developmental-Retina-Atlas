library(ggplot2)

categories <- 1:33
categories[c(11, 22, 33)] = ""

categories[1:10] = c("neurogenesis", "generation of neurons", "nervous system development", 
                "neuron differentiation", "system development", 
                "anatomical structure morphogenesis", "neuron development", 
                "multicellular organism development", 
                "cell morphogenesis involved in neuron differentiation", 
                "cell communication")

categories[12:21] = c("cytoplasmic translation", "translation", 
                  "peptide biosynthetic process", "peptide metabolic process", 
                  "amide biosynthetic process", "amide metabolic process", 
                  "organonitrogen compound biosynthetic process", 
                  "cellular nitrogen compound biosynthetic process", 
                  "protein metabolic process", "ribosome biogenesis")

categories[23:32] = c("cell adhesion", "plasma membrane bounded cell projection organization",
                " neuron development", "cell morphogenesis", " cell projection organization",
                " generation of neurons", " system development", " neurogenesis",
                " neuron differentiation", " axon guidance")

categories <- stringr::str_to_title(categories)

values <- 1:33
values[1:10] = c(3.41E-15, 7.89E-15, 8.71E-15, 3.91E-14, 2.20E-12, 
                  1.93E-11, 1.21E-10, 1.48E-10, 1.25E-09, 1.48E-09)
values[12:21] = c(8.21E-154, 3.22E-104, 8.79E-103, 2.95E-97, 1.23E-96,
                  1.34E-86, 3.11E-75, 3.57E-43, 8.73E-41, 1.81E-38)
values[23:32] = c(4.02E-10, 4.22E-10, 4.84E-10, 6.27E-10, 9.50E-10,
                  1.19E-09, 2.97E-09, 3.14E-09, 3.94E-09, 4.58E-09)

values = -log10(values)
values[c(11, 22, 33)] <- 0
group = rep(c("Module 3", "Module 2", "Module 1"), each = 11)
# Create a data frame
data <- data.frame(Category = categories, Value = values, group = group)
data <- data[33:1, ]
# Create the bar plot
ggplot(data, aes(x = Category, y = Value, fill = group)) + geom_bar(stat = "identity") +
    labs(x = "", y = "-log10(Pvalues)") + ggtitle("") +
    scale_x_discrete(limits = data$Category) + coord_flip() + theme(text = element_text(size = 30)) +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd")[3:1], limits = c("Module 3", "Module 2", "Module 1")) +
    theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent",
        color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent")) +
    theme(legend.position = c(0.8, 0.2))

ggsave("/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/PRPC_MG_Modules_GO.tiff",
    bg = "transparent", width = 20, height = 13, units = "in")