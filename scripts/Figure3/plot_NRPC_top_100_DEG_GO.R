library(ggplot2)

categories <- 1:36
categories[c(6, 12, 18, 24, 30, 36)] = ""

categories[1:5] = c("nervous system development", "    generation of neurons",
    "  neurogenesis", "    neuron differentiation", "system development")

categories[7:11] = c("regulation of postsynapse organization", "postsynapse organization",
    "cell junction organization", "regulation of synapse organization",
    "regulation of synapse structure or activity")

categories[13:17] = c("   neuron differentiation", " neuron development",
    "   generation of neurons", "neuron projection morphogenesis ", "plasma membrane bounded cell projection morphogenesis")

categories[19:23] = c(" nervous system development", " generation of neurons",
    "  neuron differentiation", " neurogenesis  ", "multicellular organism development")

categories[25:29] = c("  generation of neurons", "neurogenesis", "neuron development",
    "neuron differentiation ", "neuron projection development")

categories[31:35] = c("visual perception", "sensory perception of light stimulus",
    "multicellular organismal process", "nervous system process", "system process")
categories <- stringr::str_to_title(categories)

values <- 1:36
values[1:5] = c(4.81e-15, 1.99e-14, 3.04e-14, 2.31e-13, 3.97e-13)
values[7:11] = c(8.7e-10, 2.76e-09, 3.41e-08, 3.51e-08, 5.01e-08)
values[13:17] = c(3.91e-16, 5.88e-16, 2.43e-15, 3.02e-14, 5.42e-14)
values[19:23] = c(4.49e-11, 1.36e-09, 2.34e-09, 7.2e-09, 3.32e-08)
values[25:29] = c(2.83e-09, 2.99e-09, 9.48e-09, 2.71e-08, 1.44e-07)
values[31:35] = c(2.9e-15, 3.71e-15, 2.07e-09, 1.52e-08, 1.89e-08)

values = -log10(values)
values[c(6, 12, 18, 24, 30, 36)] <- 0
group = rep(c("RGC", "AC", "HC", "BC", "Cone", "Rod"), each = 6)
# Create a data frame
data <- data.frame(Category = categories, Value = values, group = group)
data <- data[36:1, ]
# Create the bar plot
ggplot(data, aes(x = Category, y = Value, fill = group)) + guides(fill = "none") +
    geom_bar(stat = "identity") + labs(x = "", y = "-log10(Pvalues)") +
    ggtitle("Biological Process") + scale_x_discrete(limits = data$Category) +
    coord_flip() + theme(text = element_text(size = 30)) + scale_fill_manual(values = c("#d62728",
    "#8c564b", "#2ca02c", "#bcbd22", "#e377c2", "#17becf"), limits = c("RGC",
    "AC", "HC", "BC", "Cone", "Rod")) + theme(panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent")) + theme(legend.position = c(0.8,
    0.2))

ggsave("/Users/zhenzuo/Desktop/PRPC_MG_Modules_GO.tiff", bg = "transparent",
    width = 20, height = 13, units = "in")

