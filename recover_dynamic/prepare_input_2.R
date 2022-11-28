df <- read.csv("/storage/singlecell/zz4/fetal_bash/scripts/recover_dynamic/meta.csv")

sink("/storage/singlecell/zz4/fetal_bash/scripts/recover_dynamic/recover_dynamic.sh")
for (i in 1:nrow(df)){
  cat("slurmtaco.sh --g01 -m 20G -t 1 -- python3 recover_dynamic.py ", df$INPUT[i], " ", df$OUTPUT[i], ";\n",sep = "")
}
sink()
