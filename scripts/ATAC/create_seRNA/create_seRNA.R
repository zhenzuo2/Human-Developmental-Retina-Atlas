library(ArchR)
library(Seurat)
Multi_Fetal_11w2d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_11w2d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_11w2d_FR")
Multi_Fetal_11w2d_FR_2 <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_11w2d_FR_2/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_11w2d_FR_2")
Multi_Fetal_11w2d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_11w2d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_11w2d_NR")
Multi_Fetal_13W_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_13W_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_13W_FR")
Multi_Fetal_13W_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_13W_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_13W_NR")
Multi_Fetal_14w5d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_14w5d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_14w5d_FR")
Multi_Fetal_14w5d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_14w5d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_14w5d_NR")
Multi_Fetal_19W4d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_19W4d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_19W4d_FR")
Multi_Fetal_19W4d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_19W4d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_19W4d_NR")
Multi_Fetal_20W2d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_20W2d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_20W2d_FR")
Multi_Fetal_20W2d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_20W2d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_20W2d_NR")
Multi_Fetal_23w1d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_23w1d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_23w1d_FR")
Multi_Fetal_23w1d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multi_Fetal_23w1d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multi_Fetal_23w1d_NR")
Multiome_10w_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_10w_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_10w_FR")
Multiome_10w_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_10w_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_10w_NR")
Multiome_12w3d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_12w3d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_12w3d_FR")
Multiome_12w3d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_12w3d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_12w3d_NR")
Multiome_14w2d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_14w2d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_14w2d_FR")
Multiome_14w2d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_14w2d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_14w2d_NR")
Multiome_16w4d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_16w4d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_16w4d_FR")
Multiome_16w4d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_16w4d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_16w4d_NR")
Multiome_20w1d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_20w1d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_20w1d_FR")
Multiome_20w1d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_20w1d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_20w1d_NR")
Multiome_23w4d_FR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_23w4d_FR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_23w4d_FR")
Multiome_23w4d_NR <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/Multiome_23w4d_NR/outs/filtered_feature_bc_matrix.h5",
    names = "Multiome_23w4d_NR")
sn_multiome_d59 <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/sn_multiome_d59/outs/filtered_feature_bc_matrix.h5",
    names = "sn_multiome_d59")
sn_multiome_d76c <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/sn_multiome_d76c/outs/filtered_feature_bc_matrix.h5",
    names = "sn_multiome_d76c")
sn_multiome_d76p <- import10xFeatureMatrix("/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/sn_multiome_d76p/outs/filtered_feature_bc_matrix.h5",
    names = "sn_multiome_d76p")

Multi_Fetal_11w2d_FR_2@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_11w2d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_13W_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_13W_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_14w5d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_14w5d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_19W4d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_19W4d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_20W2d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_20W2d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_23w1d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multi_Fetal_23w1d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_10w_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_10w_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_12w3d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_12w3d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_14w2d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_14w2d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_16w4d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_16w4d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_20w1d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_20w1d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_23w4d_FR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
Multiome_23w4d_NR@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
sn_multiome_d59@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
sn_multiome_d76c@rowRanges = Multi_Fetal_11w2d_FR@rowRanges
sn_multiome_d76p@rowRanges = Multi_Fetal_11w2d_FR@rowRanges

seRNA <- cbind(Multi_Fetal_11w2d_FR, Multi_Fetal_11w2d_FR_2)
seRNA <- cbind(seRNA, Multi_Fetal_11w2d_NR)
seRNA <- cbind(seRNA, Multi_Fetal_13W_FR)
seRNA <- cbind(seRNA, Multi_Fetal_13W_NR)
seRNA <- cbind(seRNA, Multi_Fetal_14w5d_FR)
seRNA <- cbind(seRNA, Multi_Fetal_14w5d_NR)
seRNA <- cbind(seRNA, Multi_Fetal_19W4d_FR)
seRNA <- cbind(seRNA, Multi_Fetal_19W4d_NR)
seRNA <- cbind(seRNA, Multi_Fetal_20W2d_FR)
seRNA <- cbind(seRNA, Multi_Fetal_20W2d_NR)
seRNA <- cbind(seRNA, Multi_Fetal_23w1d_FR)
seRNA <- cbind(seRNA, Multi_Fetal_23w1d_NR)
seRNA <- cbind(seRNA, Multiome_10w_FR)
seRNA <- cbind(seRNA, Multiome_10w_NR)
seRNA <- cbind(seRNA, Multiome_12w3d_FR)
seRNA <- cbind(seRNA, Multiome_12w3d_NR)
seRNA <- cbind(seRNA, Multiome_14w2d_FR)
seRNA <- cbind(seRNA, Multiome_14w2d_NR)
seRNA <- cbind(seRNA, Multiome_16w4d_FR)
seRNA <- cbind(seRNA, Multiome_16w4d_NR)
seRNA <- cbind(seRNA, Multiome_20w1d_FR)
seRNA <- cbind(seRNA, Multiome_20w1d_NR)
seRNA <- cbind(seRNA, Multiome_23w4d_FR)
seRNA <- cbind(seRNA, Multiome_23w4d_NR)
seRNA <- cbind(seRNA, sn_multiome_d59)
seRNA <- cbind(seRNA, sn_multiome_d76c)
seRNA <- cbind(seRNA, sn_multiome_d76p)

saveRDS(seRNA, "/storage/chentemp/zz4/adult_dev_compare/results/seRNA/seRNA.rds")
