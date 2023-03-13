df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.var.csv")
df <- df[df$fit_direction == "complete", ]

TFs <- read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_All_feature_selection_FALSE.csv")
TFs <- unique(TFs$tf)

table(df$Gene[(df$fit_direction == "complete") & (df$fit_model == 1)] %in%
        TFs)

table(df$Gene[(df$fit_direction == "complete") & (df$fit_model == 2)] %in%
        TFs)

table(df$Gene[df$fit_direction == "on"] %in% TFs)

table(df$Gene[df$fit_direction == "off"] %in% TFs)

library(ggplot2)

# create a dataset
data = data.frame(model = c(rep("Model 1", 39), rep("Model 1", 39), rep("Model 2",
                                                                        4), rep("Model 2", 15), rep("Induction", 100), rep("Induction", 23),
                            rep("Repression", 6)), Gene = c(rep("TF", 39), rep("Target gene", 39),
                                                            rep("TF", 4), rep("Target gene", 15), rep("TF", 23), rep("Target gene",
                                                                                                                     100), rep("Target gene", 6)), value = 1)
# Stacked + percent
ggplot(data, aes(fill = Gene, y = value, x = model)) + geom_bar(position = "fill",
                                                                stat = "identity") + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Percentage") + xlab("") + theme(panel.border = element_blank(),
                                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        axis.line = element_line(colour = "black")) + theme(text = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("/storage/singlecell/zz4/fetal_bash/temp/Rplot.svg")
