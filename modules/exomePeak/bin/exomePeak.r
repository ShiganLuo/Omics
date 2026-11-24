#!/usr/bin/env Rscript
# exomePeak 引擎:峰调用/差异分析 -> gene_id 翻译为 gene_name -> 清理中间物,
# 只落盘声明型输出 *_gene_names.{bed,xls}。

run_peak_calling <- function(gtf, ip_bams, input_bams, outdir) {
  exomepeak(GENE_ANNO_GTF = gtf,
            IP_BAM = ip_bams,
            INPUT_BAM = input_bams,
            OUTPUT_DIR = outdir)
}

run_diff_analysis <- function(gtf, ip_bams, input_bams, treated_ip_bams, treated_input_bams, outdir) {
  exomepeak(GENE_ANNO_GTF = gtf,
            IP_BAM = ip_bams,
            INPUT_BAM = input_bams,
            TREATED_IP_BAM = treated_ip_bams,
            TREATED_INPUT_BAM = treated_input_bams,
            OUTPUT_DIR = outdir)
}

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(".")
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
}

# 把 exomePeak 原始表 name 列的 gene_id 翻译成 gene_name 写出声明型输出,
# 全部转换成功后删除 raw_dir(原始表与 exomePeak.Rdata 均为中间物)。
annotate_and_clean <- function(raw_dir, outdir, gtf, mapping) {
  gene_env <- new.env(parent = globalenv())
  sys.source(file.path(get_script_dir(), "geneId2name.r"), envir = gene_env)
  cache_path <- file.path(outdir, "gene_map.rds")
  for (raw_name in names(mapping)) {
    infile <- file.path(raw_dir, raw_name)
    outfile <- file.path(outdir, mapping[[raw_name]])
    df <- gene_env$convert_gene_id_to_name(infile, gtf, cache_path = cache_path)
    write.table(df, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
    message("done: ", outfile)
  }
  unlink(raw_dir, recursive = TRUE)
}

main <- function() {
  suppressMessages({
    library(argparse)
    library(exomePeak)
  })
  parser <- ArgumentParser(description = "exomePeak peak calling and differential analysis with gene-name annotation")
  parser$add_argument("--gtf", required = TRUE, help = "Path to the GTF annotation file")
  parser$add_argument("--ip_bams", nargs = "+", required = TRUE, help = "List of IP BAM files")
  parser$add_argument("--input_bams", nargs = "+", required = TRUE, help = "List of Input BAM files")
  parser$add_argument("--treated_ip_bams", nargs = "*", help = "List of treated IP BAM files")
  parser$add_argument("--treated_input_bams", nargs = "*", help = "List of treated Input BAM files")
  parser$add_argument("--outdir", default = ".", help = "Output directory for results")

  args <- parser$parse_args()
  raw_dir <- file.path(args$outdir, "exomePeak_output")

  has_treated <- length(args$treated_ip_bams) > 0 && length(args$treated_input_bams) > 0
  if (has_treated) {
    run_diff_analysis(gtf = args$gtf,
                      ip_bams = args$ip_bams,
                      input_bams = args$input_bams,
                      treated_ip_bams = args$treated_ip_bams,
                      treated_input_bams = args$treated_input_bams,
                      outdir = args$outdir)
    mapping <- c("diff_peak.bed" = "diff_peaks_gene_names.bed",
                 "diff_peak.xls" = "diff_peaks_gene_names.xls",
                 "sig_diff_peak.bed" = "sig_diff_peak_gene_names.bed",
                 "sig_diff_peak.xls" = "sig_diff_peak_gene_names.xls",
                 "con_sig_diff_peak.bed" = "con_sig_diff_peak_gene_names.bed",
                 "con_sig_diff_peak.xls" = "con_sig_diff_peak_gene_names.xls")
  } else {
    run_peak_calling(gtf = args$gtf,
                     ip_bams = args$ip_bams,
                     input_bams = args$input_bams,
                     outdir = args$outdir)
    mapping <- c("peak.bed" = "all_peaks_gene_names.bed",
                 "peak.xls" = "all_peaks_gene_names.xls",
                 "con_peak.bed" = "con_peaks_gene_names.bed",
                 "con_peak.xls" = "con_peaks_gene_names.xls")
  }
  annotate_and_clean(raw_dir, args$outdir, args$gtf, mapping)
}

if (interactive() || sys.nframe() == 0) {
  main()
}
