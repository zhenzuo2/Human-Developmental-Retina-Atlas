library(ggplot2)

categories <- 1:42
categories[c(21,42)] = ""

categories[1:20] = c("synapse assembly",
                     "anterograde trans-synaptic signaling",
                     "chemical synaptic transmission",
                     "monoatomic ion transmembrane transport",
                     "trans-synaptic signaling",
                     "monoatomic ion transport",
                     "nervous system development ",
                     "synaptic signaling",
                     "cell-cell signaling",
                     "neuron development",
                     "synapse organization",
                     "modulation of chemical synaptic transmission",
                     "regulation of trans-synaptic signaling",
                     "transmembrane transport",
                     "inorganic ion transmembrane transport",
                     "aerobic electron transport chain",
                     "neuron differentiation",
                     "regulation of biological quality",
                     "cell junction organization",
                     "mitochondrial ATP synthesis coupled electron transport",
                     "ATP synthesis coupled electron transport",
                     "inorganic cation transmembrane transport")

categories[22:41] = c(
  "nervous system development",
  "neurogenesis",
  "system development",
  "generation of neurons",
  "neuron differentiation ",
  "cell junction organization ",
  "neuron development ",
  "multicellular organism development",
  "multicellular organismal process",
  "synapse organization ",
  "neuron projection development",
  "anatomical structure development",
  "postsynapse organization",
  "cell development",
  "plasma membrane bounded cell projection organization",
  "cell projection organization",
  "developmental process",
  "neuron projection morphogenesis",
  "plasma membrane bounded cell projection morphogenesis",
  "cell projection morphogenesis"
)


categories <- stringr::str_to_title(categories)

values <- 42
values[1:20] = c(1.53E-07,
                 3.26E-07,
                 3.26E-07,
                 3.28E-07,
                 3.94E-07,
                 4.56E-07,
                 6.82E-07,
                 7.00E-07,
                 2.11035E-06,
                 3.47016E-06,
                 5.24301E-06,
                 8.80653E-06,
                 9.0905E-06,
                 1.44969E-05,
                 1.82823E-05,
                 1.98363E-05,
                 2.305E-05,
                 2.66699E-05,
                 3.93691E-05,
                 4.15183E-05)
values[22:41] = c(1.67E-13,
                  4.27E-10,
                  5.09E-10,
                  4.63E-09,
                  2.87E-08,
                  2.03E-07,
                  4.55E-07,
                  4.89E-07,
                  3.38001E-06,
                  5.72407E-06,
                  1.27488E-05,
                  1.73675E-05,
                  2.07931E-05,
                  5.93226E-05,
                  7.26245E-05,
                  0.000123156,
                  0.00015163,
                  0.000167972,
                  0.000241072,
                  0.000269289)
values = -log10(values)
values[c(21, 42)] <- 0
group = rep(c("Early", "Late"), each = 21)
# Create a data frame
data <- data.frame(Category = categories, Value = values, group = group)
data = data[42:1, ]
# Create the bar plot
ggplot(data, aes(x = Category, y = Value, fill = group)) + geom_bar(stat = "identity") +
  labs(x = "", y = "-log10(Pvalues)") + ggtitle("Biological Process") +
  scale_x_discrete(limits = data$Category) + coord_flip() + theme(text = element_text(size = 30)) +
  scale_fill_manual(values = c("#1f77b4",
                               "#ff7f0e"),
                    limits = c("Early", "Late")) + theme(panel.background = element_rect(fill = "transparent"),
                                                         plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
                                                         legend.box.background = element_rect(fill = "transparent"))+ theme(legend.position = c(0.8, 0.2))

ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure3/AC_Modules_GO.svg",
       bg = "transparent", width = 20, height = 12, units = "in")