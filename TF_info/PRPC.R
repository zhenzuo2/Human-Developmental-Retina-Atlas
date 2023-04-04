TFs <- c(read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_All_feature_selection_FALSE.csv")$tf,
  read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_All_feature_selection_TRUE.csv")$tf,
  read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Macula_feature_selection_FALSE.csv")$tf,
  read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Macula_feature_selection_TRUE.csv")$tf,
  read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Peripheral_feature_selection_FALSE.csv")$tf,
  read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_Peripheral_feature_selection_TRUE.csv")$tf)

cat(unique(TFs),sep="\n")

Rod = c("ASCL1","FOXN4","HES4","HES1","HES6","NRL","OTX2")
Cone = c("ASCL1","FOXN4","HES1","HES4","OTX2","POU2F1")
BC = c("ASCL1","FEZF2","FOXN3","IRX5","IRX6","NFIA","NFIB","NFIX","OTX2","PRDM8")
AC= c("EGR1","FOXN3","FOXN4","LHX9","PRDM8")
HC= c("EGR1","FOXN4","ONECUT1","PROX1")
MG = c("GLI2","GLI3","NFIA","NFIB","NFIX","PLAGL1")
RGC = c("ASCL1","NEUROD4","EOMES","NRF1","NEUROD1","POU3F1")
