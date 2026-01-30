
load_script <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    script_dir <- if(length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
    rfiles <- list.files(
        path = file.path(script_dir, "src"),
        pattern = "\\.[Rr]$",
        full.names = TRUE,
        recursive = TRUE
    )
    log_script = file.path(script_dir, "util", "LogUtil.r")
    GeneUtil = file.path(script_dir, "util", "GeneUtil.r")
    rfiles <- c(log_script, GeneUtil, rfiles)
    lapply(rfiles, source)
    log_msg("INFO", "load", paste(rfiles, collapse = ", "))
}
load_script()
# exp.m <- read.table("/data/pub/zhousha/Totipotent20251031/RNAseqML/matrix/tpm229.tsv",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE);
log_msg("INFO", "Data loaded.")
adj.m = string_to_adj("/data/pub/zhousha/Totipotent20251031/data/STRING/9606.protein.links.v12.0.txt",
                          convert_to_entrez = TRUE,
                          species = 9606,
                          score_scale = c(0,1))
write.table(adj.m, file="/data/pub/zhousha/Totipotent20251031/data/STRING/9606_adj_matrix.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE);
log_msg("INFO", "Adjacency matrix created.")