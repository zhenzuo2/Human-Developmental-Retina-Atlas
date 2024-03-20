samples <- list.files("/storage/chentemp/zz4/adult_dev_compare/results/after_qc_seurat_object_apply_filter")

df <- as.data.frame(samples)
rownames(df) <- tools::file_path_sans_ext(df$samples)


for (x in samples) {
    print(x)
    temp <- readRDS(paste("/storage/chentemp/zz4/adult_dev_compare/results/after_qc_seurat_object_apply_filter/",
        x, sep = ""))
    x <- tools::file_path_sans_ext(x)
    df[x, "Pass Filter 1"] <- nrow(temp@meta.data)
}


for (x in samples) {
    print(x)
    temp <- readRDS(paste("/storage/chentemp/zz4/adult_dev_compare/results/DoubletFinder_seurat_object/",
        x, sep = ""))
    x <- tools::file_path_sans_ext(x)
    df[x, "Pass Filter 2"] <- nrow(temp@meta.data)
}

write.csv(df, "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table2/Number_of_Cells_Sample_part1.csv")
