args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_results_path <- args[2]

library(dplyr)
library(Seurat)
meta <- data.frame()

input_file = c("Multi_Fetal_11w2d_FR", "Multi_Fetal_11w2d_FR_2", "Multi_Fetal_11w2d_NR",
    "Multi_Fetal_13W_FR", "Multi_Fetal_13W_NR", "Multi_Fetal_14w5d_FR",
    "Multi_Fetal_14w5d_NR", "Multi_Fetal_19W4d_FR", "Multi_Fetal_19W4d_NR",
    "Multi_Fetal_20W2d_FR", "Multi_Fetal_20W2d_NR", "Multi_Fetal_23w1d_FR",
    "Multi_Fetal_23w1d_NR", "Multiome_10w_FR", "Multiome_10w_NR", "Multiome_12w3d_FR",
    "Multiome_12w3d_NR", "Multiome_14w2d_FR", "Multiome_14w2d_NR", "Multiome_16w4d_FR",
    "Multiome_16w4d_NR", "Multiome_20w1d_FR", "Multiome_20w1d_NR", "Multiome_23w4d_FR",
    "Multiome_23w4d_NR", "17W1D_Fovea_retina", "17W1D_Nasal_retina", "17W1D_Temporal_retina")

for (sam in input_file) {
    seurat_object <- readRDS(paste(input_path, sam, "_scPred.rds", sep = ""))
    meta <- bind_rows(meta, seurat_object@meta.data)
}

meta <- meta[meta$scpred_prediction %in% c("AC", "RGC", "MG", "Cone", "BC", "HC","Rod"),c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "RNA_snn_res.0.1", 
"seurat_clusters", "pANN_0.25_0.005_44", "doublet_finder", "scpred_AC", 
"scpred_Astrocyte", "scpred_BC", "scpred_Cone", "scpred_HC", 
"scpred_MG", "scpred_Microglia", "scpred_RGC", "scpred_RPE", 
"scpred_Rod", "scpred_max", "scpred_prediction", "scpred_no_rejection", 
"Time", "Region", "Days")]

meta <- meta[meta$scpred_max>0.99,]

write.table(meta, paste(output_results_path, "meta.csv", sep = ""))