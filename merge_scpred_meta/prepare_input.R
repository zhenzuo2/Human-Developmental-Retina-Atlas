input_path="/storage/singlecell/zz4/fetal_bash/results/scPred/"
df <- data.frame("Samples" = list.files(input_path,full.names = T))
df

write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/merge_scpred_meta/meta.csv",
          row.names = F)
df
