df <- read.csv("/storage/singlecell/zz4/fetal_bash/scripts/add_ldata_to_adata/meta.csv")

sink("/storage/singlecell/zz4/fetal_bash/scripts/add_ldata_to_adata/add_ldata_to_adata.sh")
for (i in 1:nrow(df)){
  cat("slurmtaco.sh -p short -m 20G -t 1 -- python3 add_ldata_to_adata.py ", df$INPUT[i], " ", df$ldata_file[i], " ", df$OUTPUT[i], ";",sep="")
  cat("\n")
}
sink()