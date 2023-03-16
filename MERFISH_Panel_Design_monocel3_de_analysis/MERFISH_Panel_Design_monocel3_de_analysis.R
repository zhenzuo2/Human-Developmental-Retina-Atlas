input_dir <- c("/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/",
    "/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/",
    "/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/",
    "/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/",
    "/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/",
    "/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/",
    "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/RPC_",
    "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_",
    "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_",
    "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/MG_")
compare_mod_region <- c(paste(input_dir, "Early_compare_mod_region.csv",
    sep = ""), paste(input_dir, "Late_compare_mod_region.csv", sep = ""),
    paste(input_dir, "All_compare_mod_region.csv", sep = ""))
fit_coefs_region <- c(paste(input_dir, "Early_fit_coefs_region_models.csv",
    sep = ""), paste(input_dir, "Late_fit_coefs_region_models.csv", sep = ""),
    paste(input_dir, "All_fit_coefs_region_models.csv", sep = ""))
res <- c()
for (i in 1:length(compare_mod_region)) {
    df <- read.csv(compare_mod_region[i])
    df <- df[complete.cases(df), ]
    df <- df[df$q_value < 0.01, ]
    df2 <- read.csv(fit_coefs_region[i])
    df2 <- df2[complete.cases(df2), ]
    df2 <- df2[(df2$q_value < 0.01) & (df2$term == "RegionPeripheral") &
        (df2$status == "OK") & (df2$normalized_effect > 1), ]
    df2 <- df2[df2$num_cells_expressed > 500, ]
    res <- c(res, intersect(df$gene_short_name, df2$gene_short_name))
}
Peripheral <- table(res)
strings <- names(table(res)[table(res) > 6])
strings <- strings[!grepl("[^[:alnum:]]", strings)]
link_strings <- grep("LINC", strings)
# remove the link strings using subsetting
if (length(link_strings) > 0) {
    strings <- strings[-link_strings]
}
# print the updated vector of strings
strings
cat(strings, sep = "\n")
deg <- data.frame(name = strings, enriched = "Peripheral")

res <- c()
for (i in 1:length(compare_mod_region)) {
    df <- read.csv(compare_mod_region[i])
    df <- df[complete.cases(df), ]
    df <- df[df$q_value < 0.01, ]
    df2 <- read.csv(fit_coefs_region[i])
    df2 <- df2[complete.cases(df2), ]
    df2 <- df2[(df2$q_value < 0.01) & (df2$term == "RegionPeripheral") &
        (df2$status == "OK") & (df2$normalized_effect < -1), ]
    df2 <- df2[df2$num_cells_expressed > 500, ]
    res <- c(res, intersect(df$gene_short_name, df2$gene_short_name))
}
Macula <- table(res)
strings <- names(table(res)[table(res) > 5])
strings <- strings[!grepl("[^[:alnum:]]", strings)]
if (length(link_strings) > 0) {
    strings <- strings[-link_strings]
}
# print the updated vector of strings
strings

deg <- rbind(deg, data.frame(name = strings, enriched = "Macula"))

genes = c("RARA", "RARB", "RARG", "CYP26A1", "CYP26B1", "CYP26C1", "CRABP1",
    "CRABP2", "DHRS3", "BMP2", "BMP4", "BMP7", "FGF1", "FGF2", "FGF3",
    "FGF8", "PTCH1", "ALDH1A1", "CRMP1", "FGFR4", "VEGF", "GFAP", "KIT",
    "CD44", "COL4A3")
for (gene in genes) {
    if (!(gene %in% deg$name)) {
        print(gene)
        print(Peripheral[gene])
        print(Macula[gene])
    }
}
deg <- rbind(deg, c("RARB", "Macula"))
deg <- rbind(deg, c("FGF1", "Macula"))
deg <- rbind(deg, c("CRABP1", "Macula"))
deg <- rbind(deg, c("CRABP2", "Peripheral"))
deg <- rbind(deg, c("DHRS3", "Peripheral"))
deg <- rbind(deg, c("BMP7", "Macula"))

deg <- deg[order(deg$enriched, deg$name), ]

cat(deg$name, sep = "\n")

labs <- tools::file_path_sans_ext(paste(sapply(strsplit(compare_mod_region,
    "/"), "[[", 7), sapply(strsplit(compare_mod_region, "/"), "[[", 8),
    sep = "_"))
deg[, labs] <- FALSE
for (i in 1:length(compare_mod_region)) {
    df <- read.csv(compare_mod_region[i])
    df <- df[complete.cases(df), ]
    df <- df[df$q_value < 0.01, ]
    deg[, labs[i]] <- deg$name %in% df$gene_id
    print(labs[i])
    print(nrow(df))
}

labs <- tools::file_path_sans_ext(paste(sapply(strsplit(fit_coefs_region,
    "/"), "[[", 7), sapply(strsplit(fit_coefs_region, "/"), "[[", 8), sep = "_"))
deg[, labs] <- FALSE
for (i in 1:length(fit_coefs_region)) {
    df2 <- read.csv(fit_coefs_region[i])
    df2 <- df2[complete.cases(df2), ]
    df2 <- df2[(df2$q_value < 0.01) & (df2$term == "RegionPeripheral") &
        (df2$status == "OK") & (df2$normalized_effect < -1), ]
    df2 <- df2[df2$num_cells_expressed > 100, ]
    deg[, labs[i]] <- deg$name %in% df$gene_id
    print(labs[i])
    print(nrow(df2))
}
######################################################################################################################################################

input_dir <- c("/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/RPC_",
    "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_",
    "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_")
compare_mod_region <- c(paste(input_dir, "Early_compare_mod_region.csv",
    sep = ""), paste(input_dir, "Late_compare_mod_region.csv", sep = ""),
    paste(input_dir, "All_compare_mod_region.csv", sep = ""))
fit_coefs_region <- c(paste(input_dir, "Early_fit_coefs_region_models.csv",
    sep = ""), paste(input_dir, "Late_fit_coefs_region_models.csv", sep = ""),
    paste(input_dir, "All_fit_coefs_region_models.csv", sep = ""))
res <- c()
for (i in 1:length(compare_mod_region)) {
    df <- read.csv(compare_mod_region[i])
    df <- df[complete.cases(df), ]
    df <- df[df$q_value < 0.01, ]
    df2 <- read.csv(fit_coefs_region[i])
    df2 <- df2[complete.cases(df2), ]
    df2 <- df2[(df2$q_value < 0.01) & (df2$term == "RegionPeripheral") &
        (df2$status == "OK") & (df2$normalized_effect > 1), ]
    df2 <- df2[df2$num_cells_expressed > 100, ]
    res <- c(res, intersect(df$gene_short_name, df2$gene_short_name))
}
Peripheral <- table(res)
strings <- names(table(res)[table(res) > 5])
strings <- strings[!grepl("[^[:alnum:]]", strings)]
link_strings <- grep("LINC", strings)
# remove the link strings using subsetting
if (length(link_strings) > 0) {
    strings <- strings[-link_strings]
}
# print the updated vector of strings
strings
cat(strings, sep = "\n")
deg2 <- data.frame(name = strings, enriched = "Peripheral")

res <- c()
for (i in 1:length(compare_mod_region)) {
    df <- read.csv(compare_mod_region[i])
    df <- df[complete.cases(df), ]
    df <- df[df$q_value < 0.01, ]
    df2 <- read.csv(fit_coefs_region[i])
    df2 <- df2[complete.cases(df2), ]
    df2 <- df2[(df2$q_value < 0.01) & (df2$term == "RegionPeripheral") &
        (df2$status == "OK") & (df2$normalized_effect < -1), ]
    df2 <- df2[df2$num_cells_expressed > 100, ]
    res <- c(res, intersect(df$gene_short_name, df2$gene_short_name))
}
Macula <- table(res)
strings <- names(table(res)[table(res) > 2])
strings <- strings[!grepl("[^[:alnum:]]", strings)]
link_strings <- grep("LINC", strings)
# remove the link strings using subsetting
if (length(link_strings) > 0) {
    strings <- strings[-link_strings]
}
# print the updated vector of strings
strings

deg2 <- rbind(deg2, data.frame(name = strings, enriched = "Macula"))

deg2 <- deg[order(deg2$enriched, deg$name), ]

cat(deg2$name, sep = "\n")

table(deg2$enriched)

dput(deg[deg$enriched == "Macula", "name"])

dput(deg[deg$enriched == "Peripheral", "name"])