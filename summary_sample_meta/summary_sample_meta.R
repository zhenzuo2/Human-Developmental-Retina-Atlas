mapping <- read.csv("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal_sample_meta.csv")
samples <- list.dirs("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                     full.names = F, recursive = F)
df1 = data.frame()
df2 = data.frame()
Sample.ID = c()
for (sample in samples) {
  print(sample)
  if (file.exists(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                        sample, "/outs/metrics_summary.csv", sep = ""))) {
    df1 <- plyr::rbind.fill(df1, read.csv(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                                              sample, "/outs/metrics_summary.csv", sep = "")))
    Sample.ID = c(Sample.ID,sample)
    df1$Sample.ID = Sample.ID
  } else {
    df2 <- plyr::rbind.fill(df2, read.csv(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                                              sample, "/outs/summary.csv", sep = "")))
  }
}


write.table(df1, "/storage/singlecell/zz4/fetal_bash/results/meta/snRNA-seq_meta.csv",
            quote = T, sep = ",", row.names = FALSE)

write.table(df2, "/storage/singlecell/zz4/fetal_bash/results/meta/multiomics_meta.csv",
            quote = T, sep = ",", row.names = FALSE)


adata.obs["Weeks"] = adata.obs.Days.map(
    {
        70: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW16",
        103: "PCW16",
        116: "PCW16",
        120: "PCW16",
        136: "PCW20",
        137: "PCW20",
        141: "PCW20",
        142: "PCW20",
        162: "PCW23",
        165: "PCW23",
    }
)