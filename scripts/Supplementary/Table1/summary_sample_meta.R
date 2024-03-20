
mapping <- read.csv("/storage/chentemp/zz4/adult_dev_compare/data/Sample_meta/Retina_fetal_sample_meta.csv")
samples <- list.dirs("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/",
    full.names = F, recursive = F)
df = data.frame()
Sample.ID = c()
for (sample in samples) {
    print(sample)
    df <- plyr::rbind.fill(df, read.csv(paste("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/",
        sample, "/outs/summary.csv", sep = "")))
}
df$Sample.ID <- samples
df$Time <- mapping$Time
df$Region <- mapping$Region
df$Days <- mapping$Days
df$Data.Type <- mapping$Data.Type

from = c(59, 70, 76, 79, 87, 91, 100, 103, 116, 137, 141, 142, 162, 165)

to = c("PCW8", "PCW10", "PCW10", "PCW10", "PCW13", "PCW13", "PCW15", "PCW15",
    "PCW15", "PCW19", "PCW19", "PCW19", "PCW23", "PCW23")
df$Weeks.Group <- plyr::mapvalues(df$Days, from = from, to = to)

df = df[c(1, 49:53, 2:48)]

write.table(df, "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table1/multiomics_meta.csv",
    quote = T, sep = ",", row.names = FALSE)
