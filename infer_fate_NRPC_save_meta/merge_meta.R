
AC<- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_AC.csv")
HC<- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_HC.csv")
RGC <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_RGC.csv")
BC<- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_BC.csv")
Rod<- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Rod.csv")
Cone<- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Cone.csv")

#put all data frames into list
df_list <- list(AC, HC, RGC)
#merge all data frames in list
df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_AC_HC_RGC.csv")

#put all data frames into list
df_list <- list(BC, Rod, Cone)
#merge all data frames in list
df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_BC_Rod_Cone.csv")
