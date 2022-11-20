df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/RPC_major_sub_class_obs.csv")
df <-df[df$subclass=="NRPC",]
colnames(df)[1]=""
write.csv(df, "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/NRPC_major_sub_class_obs.csv",row.names = F)

df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/RPC_major_sub_class_obs.csv")
df <-df[df$subclass=="PRPC",]
colnames(df)[1]=""
write.csv(df,"/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/PRPC_major_sub_class_obs.csv",row.names = F)

df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/RPC_major_sub_class_obs.csv")
df <-df[df$subclass=="MG",]
colnames(df)[1]=""
write.csv(df,"/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/MG_major_sub_class_obs.csv",row.names = F)

input_file_path=c(
  "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label/AC_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label/BC_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label/Cone_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label/Rod_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label/HC_major_sub_class_obs.csv",
    "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label/RGC_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/RPC_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/NRPC_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/PRPC_major_sub_class_obs.csv",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/MG_major_sub_class_obs.csv"
)


output_file_path=c(
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_AC/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_BC/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_Cone/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_Rod/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_HC/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_RGC/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_RPC/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_NRPC/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_PRPC/",
  "/storage/singlecell/zz4/fetal_bash/results/pseudotime_MG/"
)
df <- data.frame("INPUT" = input_file_path, "OUTPUT"=output_file_path)
write.table(df,"/storage/singlecell/zz4/fetal_bash/scripts/monocle3_pseudotime/meta.csv",row.names = F)
