pseudotime_file <- "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv"
time <- read.csv(pseudotime_file)

time$Days_ <- mapvalues(time$Days, from = c(70, 79, 87, 91, 100, 103, 116,
    136, 137, 141, 142, 162, 165), to = c("FW10", "FW10", "FW13", "FW13",
    "FW13", "FW16", "FW16", "FW19", "FW19", "FW19", "FW19", "FW23", "FW23"),
    warn_missing = TRUE)


library(ggplot2)
# grouped boxplot
p <- ggplot(time, aes(x = Days_, y = latent_time, fill = Region)) + geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = quantile(data$y, c(0.1, 0.9))) + xlab("Clock time") +
    ylab("Latent time") + theme_bw() + theme(text = element_text(size = 21))
svg("/storage/singlecell/zz4/fetal_bash/figures/Latent_time_vs_clock_time/PRPC.svg")
p
dev.off()