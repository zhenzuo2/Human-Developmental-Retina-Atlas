df1 <- read.csv("/storage/chentemp/zz4/adult_dev_compare/data/Sample_meta/Retina_fetal_sample_meta.csv")
df2 <- read.csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table2/Number_of_Cells_Sample_part2.csv")
rownames(df1) <- df1$Samples
rownames(df2) <- df2$smapleID
df <- merge(df1, df2, by.x = 0, by.y = 0)
colnames(df)
df <- df[, c("Samples", "Days", "Number.of.Cells.Sequenced", "Pass.Filter.1",
    "Pass.Filter.2", "Pass.Filter.3", "Pass.Filter.4")]
df <- df[order(df$Days), ]
colnames(df)[3] <- "Number.of.Nuclei.Sequenced"

write.table(df, "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table2/Number_of_Cells_Sample_part3.csv",
    quote = F, sep = ",", row.names = F)