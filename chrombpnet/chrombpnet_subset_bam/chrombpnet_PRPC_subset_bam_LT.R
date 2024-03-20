subsetbam <- function(bam, cell_id, name, sh_folder, output) {
    sink(paste(sh_folder, name, ".sh", sep = ""))
    cat(paste("echo \"Processing ", name, "\"\n", sep = ""))
    cat(paste("samtools view -H ", bam, " > ", output, name, "_SAM_header",
        sep = ""))
    cat("\n")
    cat(paste("samtools view ", bam, " | LC_ALL=C grep -F -f ", cell_id,
        " > ", output, name, sep = ""))
    cat("\n")
    cat(paste("cat ", output, name, "_SAM_header ", output, name, " > ",
        output, name, ".sam", sep = ""))
    cat("\n")
    cat(paste("samtools view -b ", output, name, ".sam > ", output, name,
        ".bam", sep = ""))
    cat("\n")
    sink()
}

for (label in c("PRPC_1", "PRPC_2", "PRPC_3", "PRPC_4")) {
    sh_folder = paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_subsetbam/",
        label, "/", sep = "")
    dir.create(sh_folder, showWarnings = FALSE)
    for (f in list.files(paste("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/",
        label, sep = ""), full.names = T)) {
        print(f)
        name = print(tools::file_path_sans_ext(basename(f)))
        bam = paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
            name, "/outs/atac_possorted_bam.bam", sep = "")
        cell_id = f
        output = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_subsetbam/",
            label, "/", sep = "")
        dir.create(output, showWarnings = FALSE)
        subsetbam(bam, cell_id, name, sh_folder, output)
    }
}

## cd /storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_subsetbam/PRPC_1
## for f in *.sh; do slurmtaco.sh --g00 -m 30G -t 1 -- sh $f; done
## for f in *.sh; do slurmtaco.sh --d01 -m 10G -t 1 -- sh $f; done
